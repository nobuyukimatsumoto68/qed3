#!/usr/bin/env python3
# hankel_rebase_scan_claude.py
#   Scan the block-Hankel + rebase knobs of the 0++ two-sigma / two-meson GEVP.
#   Reuses per_config / build_matrix_conn / gevp / rebase_vectors / rebased_effmass from the
#   working driver gevp_twosigma_interacting_claude.py (imported as drv, NOT modified).
#
#   Economy: the expensive per-config correlator store is computed ONCE and cached to the
#   scratchpad as a pickle; every knob combination is then evaluated in-process (cheap linear
#   algebra), so the whole scan costs one store-build (~50 s) instead of one per configuration.
#
#   NEW code paths added here (do not exist in the driver):
#     - hankel_off(Cts, offsets): generalized block-Hankel with an EXPLICIT offset list
#       Chat(t)_{(a,i),(b,j)} = C(t + off[a] + off[b]); supports non-uniform ladders
#       (Dt=1,2,3 -> offsets 0,1,2,3 ; Dt=1..5 -> offsets 0,1,2,3,4,5). t_max = twin - 2*max(off).
#     - staged (multiple) rebase: compose per-stage rebase projections V1 @ V2 @ ...
#       (algebraically identical to successive project-then-rebase, see note in staged_project).
#
#   Fixed ensemble:
#     ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
#     NVDIR=distill_Nv24 , SPLIT=1 , CONN=1 , BINSIZE=10 (jackknife blocking for autocorrelation).
#
#   Refs: block-Hankel / GPOF Aubin-Orginos 1010.0202 ; distillation Peardon 0905.2160 ;
#         GEVP Z-factors Blossier 0902.1265 ; sigma^2-F^2 mixing Chester-Pufu 1603.05582.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import pickle
import numpy as np
import distill_contract_claude as dc
import gevp_twosigma_interacting_claude as drv

SPLIT = 1
BINSIZE = 10
SCRATCH = "/tmp/claude-1000/-mnt-barracuda22-qed3/cdb76ded-4ed0-4f78-8eda-57542aa97f37/scratchpad"


def build_store():
    # compute the per-config correlator store once (mirrors drv.main's loop)
    tag = dc.ENS.split("nu0")[0]
    KS = dc.KS
    ncfg = len(KS)
    cache = "%s/store_%s_%d_claude.pkl" % (SCRATCH, tag.replace(".", "p"), ncfg)
    if os.path.exists(cache):
        with open(cache, "rb") as f:
            store, twin = pickle.load(f)
        print("# loaded cached store  ncfg=%d  twin=%d  <- %s" % (ncfg, twin, cache))
        return store, twin, ncfg, tag
    print("# building store (one-time)  ncfg=%d ..." % ncfg)
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    keys = ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]
    store = {k: [] for k in keys}
    twin = None
    for j, k in enumerate(KS):
        diags, CAA, Css, C2s, C2A, o2, oA, o2s, tw = drv.per_config(k, w00, SPLIT)
        twin = tw
        store["diags"].append(diags)
        store["CAA"].append(CAA)
        store["Css"].append(Css)
        store["C2s"].append(C2s)
        store["C2A"].append(C2A)
        store["o2"].append(o2)
        store["oA"].append(oA)
        store["o2s"].append(o2s)
        if (j + 1) % 100 == 0:
            print("#   ... %d/%d" % (j + 1, ncfg))
    with open(cache, "wb") as f:
        pickle.dump((store, twin), f)
    print("# cached store -> %s" % cache)
    return store, twin, ncfg, tag


# ---------------------------------------------------------------------------
# generalized block-Hankel with an explicit offset list (NEW)
# ---------------------------------------------------------------------------
def hankel_off(Cts, offsets):
    twin, N, _ = Cts.shape
    nb = len(offsets)
    omax = max(offsets)
    tmax = twin - 2 * omax
    Big = np.full((tmax, nb * N, nb * N), np.nan)
    for t in range(tmax):
        for a in range(nb):
            for b in range(nb):
                Big[t, a * N:(a + 1) * N, b * N:(b + 1) * N] = Cts[t + offsets[a] + offsets[b]]
    return Big


def hankel_permode(Cts, offs_list):
    # PER-OPERATOR block-Hankel: operator i carries its own offset list offs_list[i].
    # Augmented dim D = sum_i len(offs_list[i]); rows/cols enumerate (i, p) for p in offs_list[i].
    # Block ((i,p),(j,q)) = Cts[t + p + q][i, j].  (hankel_off is the special case of one common list.)
    twin, N, _ = Cts.shape
    assert len(offs_list) == N, "offs_list must have one entry per operator"
    rows = []
    for i in range(N):
        for p in offs_list[i]:
            rows.append((i, p))
    D = len(rows)
    omax = max(max(o) for o in offs_list)
    tmax = twin - 2 * omax
    Big = np.full((tmax, D, D), np.nan)
    for t in range(tmax):
        for r, (i, p) in enumerate(rows):
            for c, (j, q) in enumerate(rows):
                Big[t, r, c] = Cts[t + p + q][i, j]
    return Big


def staged_project(Big, stages, t0):
    # stages = list of (t_reb, nkeep) applied successively.  Returns the composite projection
    # Vtot (Nbig x final_nkeep).  Composition is exact: with Cr1 = V1^T Big V1 and
    # V2 = rebase(Cr1), Cr2 = V2^T Cr1 V2 = (V1 V2)^T Big (V1 V2), so Vtot = V1 @ V2 @ ...
    Vtot = None
    cur = Big
    for (tr, nk) in stages:
        V = drv.rebase_vectors(cur, tr, t0, nk)     # (dim_cur x nk)
        if Vtot is None:
            Vtot = V
        else:
            Vtot = Vtot @ V
        cur = np.einsum("ai,tab,bj->tij", V, cur, V)
    return Vtot


def rebased_effmass_fixed(Big, Vop, t0):
    # like drv.rebased_effmass but ALWAYS returns (tmax-1, nkeep): levels the reduced metric
    # cannot resolve (rank-deficient jackknife samples, e.g. a fine ladder on a rank-2 base) are
    # padded with nan so the jackknife arrays stay homogeneous.
    Cr = np.einsum("ai,tab,bj->tij", Vop, Big, Vop)
    twinr, nk, _ = Cr.shape
    C0 = 0.5 * (Cr[t0] + Cr[t0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > 1e-10 * wv.max()
    nlev = int(keep.sum())
    Uk = Uv[:, keep] / np.sqrt(np.abs(wv[keep]))
    lam = np.full((twinr, nk), np.nan)
    for t in range(twinr):
        Mt = Uk.T @ (0.5 * (Cr[t] + Cr[t].T)) @ Uk
        try:
            lam[t, :nlev] = np.sort(np.linalg.eigvals(Mt).real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        return np.log(lam[:-1] / lam[1:])


def jackknife_bin(store, keys, binsize):
    ncfg = len(store["o2"])
    nbin = ncfg // binsize
    jkstore = {key: [np.mean([store[key][i * binsize + r] for r in range(binsize)], axis=0)
                     for i in range(nbin)] for key in keys}
    return jkstore, nbin


KEYS = ["diags", "CAA", "Css", "C2s", "C2A", "o2", "oA", "o2s"]


def eval_config(store, twin, nooa, offsets, stages, t0):
    # build central + jackknife Hankel matrices, apply the (staged) rebase, return effmass table
    drv.NOOA = nooa
    ncfg = len(store["o2"])
    Dcen = drv.ens_avg(store, list(range(ncfg)))
    Ccen = drv.build_matrix_conn(Dcen, twin)
    Big = hankel_off(Ccen, offsets)
    Vtot = staged_project(Big, stages, t0)
    em_c = rebased_effmass_fixed(Big, Vtot, t0)        # (tmax-1, nkeep)
    jkstore, nbin = jackknife_bin(store, KEYS, BINSIZE)
    ems = []
    for i in range(nbin):
        idx = [j for j in range(nbin) if j != i]
        Bi = hankel_off(drv.build_matrix_conn(drv.ens_avg(jkstore, idx), twin), offsets)
        ems.append(rebased_effmass_fixed(Bi, Vtot, t0))  # FIXED central Vtot
    ems = np.array(ems)                                # (nbin, tmax-1, nkeep)
    em_err = np.sqrt((nbin - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    return em_c, em_err, ems, nbin


def window_mass(em_c, ems, nbin, lev, tlo, thi):
    # jackknife window-average mass + error for one level over [tlo, thi] (inclusive t-indices)
    sl = slice(tlo, thi + 1)
    m_c = np.nanmean(em_c[sl, lev])
    mj = np.nanmean(ems[:, sl, lev], axis=1)
    e = np.sqrt((nbin - 1) * np.mean((mj - mj.mean()) ** 2))
    return m_c, e


def window_slope(em_c, ems, nbin, lev, tlo, thi):
    # jackknife slope (per unit t) of the effmass over [tlo, thi] via least squares
    ts = np.arange(tlo, thi + 1)
    ok = np.isfinite(em_c[tlo:thi + 1, lev])
    ts = ts[ok]
    if len(ts) < 3:
        return np.nan, np.nan
    A = np.vstack([ts, np.ones_like(ts)]).T
    y = em_c[tlo:thi + 1, lev][ok]
    s_c = np.linalg.lstsq(A, y, rcond=None)[0][0]
    sj = []
    for i in range(nbin):
        yi = ems[i, tlo:thi + 1, lev][ok]
        sj.append(np.linalg.lstsq(A, yi, rcond=None)[0][0])
    sj = np.array(sj)
    e = np.sqrt((nbin - 1) * np.mean((sj - sj.mean()) ** 2))
    return s_c, e


if __name__ == "__main__":
    store, twin, ncfg, tag = build_store()
    print("# store ready: twin=%d ncfg=%d tag=%s" % (twin, ncfg, tag))
