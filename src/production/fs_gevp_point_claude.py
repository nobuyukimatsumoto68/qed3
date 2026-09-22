#!/usr/bin/env python3
# fs_gevp_point_claude.py  [CHUNK 2/3: point-operator connected GEVP {sigma^2_00, O_2m, O_1m}, FS]
# Run:  ENS=... NVDIR=distill_Nv24 SPLIT=1 T0=3 BINSIZE=10 VALIDATE=1 python3 fs_gevp_point_claude.py
#
# Operators (all sigma=sigma_FS density; geometry projects onto the wanted state; full-connected A-G):
#   0 sigma^2_00 : K[i,j]=(A_i Y00)(A_j Y00) , times (T,T)      -- Y00 primary target
#   1 O_2m       : K[i,j]=A_i delta_{j,P(i)} , times (T,T)      -- antipodal -> two-meson
#   2 O_1m       : K[i,j]=A_i delta_{i,j}    , times (T,T+delta)-- coincident split -> one-meson
# Connected correlator <O_i(T) O_j(S)> = sum over BRIDGING perms (sink{0,1}<->source{2,3}) of
#   (-1)^#cyc * einsum[ prod cycle A-blocks (spin-traced) , Ksink[i,j], Ksrc[k,l] ] , summed leg S(tau)+Stilde(-tau').
# Position propagator block A_leg(ta,tb) = U_ta leg(ta,tb) U_tb^dag reshaped (nsite,2,nsite,2);
#   S : A(ta,tb), diagonal ta==tb -> A - 1/2 I (ultralocal GW contact).
#   Stilde: -A'(ta,tb) [A'=U tau' U^dag] ; ta==tb -> -1/2 (A(t,t)+A'(t,t)) [furnished contact 1/2(tau'-tau)].
# VALIDATE=1: C[0][0] (Y00) connected must equal diag_effmass diags_pair connected (A,B,C,D,E,G) sum.
# See fs_gevp_connected_impl_plan_claude.md (chunks).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
from itertools import permutations
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

NS = dc.NS
SPLIT = int(os.environ.get("SPLIT", "1"))
VALIDATE = int(os.environ.get("VALIDATE", "1"))
# MODE_CONTACT (env, default 0): where the GW equal-time contact 1/2 is subtracted.
#   0 (original) = -1/2 I on the FULL 2Ns space (A -= 0.5 Iv AFTER U tau U^dag).  EXACT only at the COMPLETE
#     basis (V V^dag = I); under TRUNCATION (V V^dag = P != I) the contact on the removed modes is uncompensated,
#     which REVIVES the single-sigma tadpole D_S = Tr[Phi(tau_tt - 1/2)] (grows linearly in # removed modes) and
#     feeds a single-meson piece into <sigma^2 sigma^2>, collapsing C_22 toward m_PS.
#   1 (fix)      = -1/2 in MODE space: subtract 1/2 I_Nv from the Nv x Nv peram BEFORE U(...)U^dag (= -1/2 P in
#     position space).  The GW contact is diagonal (=1/2) in the distillation basis, so this is preserved exactly
#     under slicing -> D_S stays 0 truncated (verified tadpole_trunc_check_claude.py).  Identical to 0 at complete.
MODE_CONTACT = int(os.environ.get("MODE_CONTACT", "0"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
DTMAX = int(os.environ.get("DTMAX", "24"))
NCFG = int(os.environ.get("NCFG", "0"))          # 0 = all
TIMING = int(os.environ.get("TIMING", "0"))
NPROC = int(os.environ.get("NPROC", "1"))        # config-parallel workers (each OMP=1)

# ============================ MACRO: FS "Stilde" furnished leg (DISABLED) ============================
# FS_STILDE_FURNISH_ENABLED = False (default) GUARDS OFF every CALL SITE of the FS Stilde furnished leg
#   AblkSt -- so AblkSt is never invoked (the guard lives at the call sites, not inside AblkSt).  At m=0 the
#   Stilde leg collapses to the plain S leg (tau, = AblkS) by the Ginsparg-Wilson identity, so sigma_FS ==
#   sigma_PS and callers use the S-part only (the S+Stilde sum becomes 2 x the S-part).
#   The old AblkSt built the stored tau_gw = -(1-D_ov^dag)D_ov^{-1}, which applies (1-D_ov^dag) to the
#   FORWARD inverse and carries a spurious BARE D_ov^dag -- an O(1) ARTIFACT (it faked an FS single-meson
#   coupling).  Full derivation + why PS==FS: fs_gw_collapse_v_agent_two_meson_claude.md (this dir).
# Flip to True ONLY for a CORRECT massive (m != 0) implementation, where GW is modified, the collapse
#   fails, and a GENUINE Stilde leg is needed -- but that leg must use the ADJOINT inverse D_ov^{-dag}
#   (= delta - tau), NOT the old forward-furnished tau_gw body.  See AblkSt below.
FS_STILDE_FURNISH_ENABLED = False
# ====================================================================================================

# bridging permutations (sink {0,1} <-> source {2,3}) with their cycle decomposition
PERMS = []
for pi in permutations(range(4)):
    seen = [False] * 4
    cycles = []
    for i in range(4):
        if seen[i]:
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(j)
            j = pi[j]
        cycles.append(cyc)
    conn = any(any(v in (0, 1) for v in c) and any(v in (2, 3) for v in c) for c in cycles)
    if conn:
        PERMS.append(cycles)

SITE = ['i', 'j', 'k', 'l']
SPIN = ['A', 'B', 'C', 'D']


def antipodal_map():
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    n = sites.shape[0]
    P = np.array([int(np.argmin(np.linalg.norm(sites + sites[i], axis=1))) for i in range(n)])
    return P


_PATH_CACHE = {}


def perm_contrib(cycles, tvec, Ablk, Ksink, Ksrc):
    arrs = []
    idxs = []
    for cyc in cycles:
        m = len(cyc)
        for a in range(m):
            va = cyc[a]
            vb = cyc[(a + 1) % m]
            arrs.append(Ablk(tvec[va], tvec[vb]))
            idxs.append(SITE[va] + SPIN[va] + SITE[vb] + SPIN[vb])
    arrs.append(Ksink)
    idxs.append('ij')
    arrs.append(Ksrc)
    idxs.append('kl')
    sub = ','.join(idxs) + '->'
    path = _PATH_CACHE.get(sub)
    if path is None:
        path = np.einsum_path(sub, *arrs, optimize='optimal')[0]
        _PATH_CACHE[sub] = path
    val = np.einsum(sub, *arrs, optimize=path)
    return ((-1.0) ** len(cycles)) * val


def make_config(k):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    n2 = V.shape[2] if V.ndim == 3 else None
    U = [V[tsrc0 + a].T for a in range(twin)]     # (2Ns, Nv)
    twoNs = U[0].shape[0]
    nsite = twoNs // NS
    Iv = np.eye(twoNs)
    cacheS = {}
    cacheSt = {}

    def AblkS(ta, tb):
        key = (ta, tb)
        if key not in cacheS:
            if ta == tb and MODE_CONTACT:
                A = U[ta] @ (tau[ta, tb] - 0.5 * np.eye(tau.shape[-1])) @ U[tb].conj().T
            else:
                A = U[ta] @ tau[ta, tb] @ U[tb].conj().T
                if ta == tb:
                    A = A - 0.5 * Iv
            cacheS[key] = A.reshape(nsite, NS, nsite, NS)
        return cacheS[key]

    # AblkSt = the FS "Stilde" furnished leg.  It is DISABLED AT THE CALL SITES via the
    #   FS_STILDE_FURNISH_ENABLED macro (near the top of this file; default False), so this function is NEVER
    #   INVOKED in the default (massless) program.  Reason: at m=0 the GW collapse makes this leg equal to the
    #   plain S leg AblkS (S~ D_ov^{-dag} = 1 - D_ov^{-dag} = D_ov^{-1} = tau), so sigma_FS == sigma_PS and the
    #   callers use AblkS only (the S+Stilde sum becomes 2 x the S-part).  The body below is the OLD furnished
    #   form -(1-D_ov^dag)D_ov^{-1} = the stored tau_gw, which applies (1-D_ov^dag) to the FORWARD inverse and
    #   carries a spurious bare D_ov^dag = an O(1) artifact (it faked an FS single-meson coupling).  Full
    #   derivation + why PS==FS: fs_gw_collapse_v_agent_two_meson_claude.md (this dir).  Flipping the macro to
    #   True re-invokes THIS body; for massive m != 0 (the only case a genuine Stilde leg is needed) it must
    #   FIRST be rewritten to use the ADJOINT inverse D_ov^{-dag} = delta - tau, NOT this forward-furnished form.
    # >>> COMMAND: DO NOT DELETE THIS COMMENT OR THE AblkSt BODY.  They record a fixed bug (the tau_gw FS-leg
    #     artifact, disabled at the call sites) and the massive-case TODO (a correct adjoint delta-tau leg). <<<
    def AblkSt(ta, tb):
        key = (ta, tb)
        if key not in cacheSt:
            if ta == tb:
                A = -0.5 * (U[ta] @ (taugw[ta, ta] + tau[ta, ta]) @ U[ta].conj().T)
            else:
                A = -(U[ta] @ taugw[ta, tb] @ U[tb].conj().T)
            cacheSt[key] = A.reshape(nsite, NS, nsite, NS)
        return cacheSt[key]

    return AblkS, AblkSt, twin, nsite, U, tau


def make_config_win(k, w):
    # Sibling of make_config for a SPECIFIC source window w (nsrc>1 data).  Builds AblkS/AblkSt from that
    # window's own tau_w / taugw_w and its base source timeslice tsrc0_w (V is shared all-t).  NO cross-window
    # object is ever formed.  make_config_win(k, 0) on any file is IDENTICAL to make_config(k) (window 0).
    # Returns the same tuple as make_config plus tsrc0_w (so callers index global objects, e.g. F^2, per window).
    V, windows = dc.load_peram_windows(k)
    tsrc0, tau, taugw = windows[w]
    twin = tau.shape[0]
    U = [V[tsrc0 + a].T for a in range(twin)]     # (2Ns, Nv), indexed at this window's absolute times
    twoNs = U[0].shape[0]
    nsite = twoNs // NS
    Iv = np.eye(twoNs)
    cacheS = {}
    cacheSt = {}

    def AblkS(ta, tb):
        key = (ta, tb)
        if key not in cacheS:
            if ta == tb and MODE_CONTACT:
                A = U[ta] @ (tau[ta, tb] - 0.5 * np.eye(tau.shape[-1])) @ U[tb].conj().T
            else:
                A = U[ta] @ tau[ta, tb] @ U[tb].conj().T
                if ta == tb:
                    A = A - 0.5 * Iv
            cacheS[key] = A.reshape(nsite, NS, nsite, NS)
        return cacheS[key]

    # AblkSt = the FS "Stilde" furnished leg (window sibling of make_config's).  DISABLED AT THE CALL SITES
    #   via the FS_STILDE_FURNISH_ENABLED macro (default False) -- never invoked in the massless program; at
    #   m=0 it collapses to AblkS (tau) by GW, so sigma_FS == sigma_PS.  See make_config's AblkSt above and
    #   fs_gw_collapse_v_agent_two_meson_claude.md for the full rationale (buggy forward-furnished tau_gw =
    #   spurious bare D_ov^dag; massive m != 0 needs the ADJOINT delta-tau leg, not this body).
    # >>> COMMAND: DO NOT DELETE THIS COMMENT OR THE AblkSt BODY. <<<
    def AblkSt(ta, tb):
        key = (ta, tb)
        if key not in cacheSt:
            if ta == tb:
                A = -0.5 * (U[ta] @ (taugw[ta, ta] + tau[ta, ta]) @ U[ta].conj().T)
            else:
                A = -(U[ta] @ taugw[ta, tb] @ U[tb].conj().T)
            cacheSt[key] = A.reshape(nsite, NS, nsite, NS)
        return cacheSt[key]

    return AblkS, AblkSt, twin, nsite, U, tau, tsrc0


def op_vspec(op, letters, dual, w):
    # per-vertex (site-letter, weight-vector-or-None, antipode-reindex) for the operator's 2 vertices
    l0, l1 = letters
    if op == 0:                         # sigma^2_00 : separable w_i w_j  (two free indices, each weighted w)
        return [(l0, w, False), (l1, w, False)]
    if op == 1:                         # O_2m antipodal : v1 site = P(v0) (shared index, antipode-reindex)
        return [(l0, dual, False), (l0, None, True)]
    return [(l0, dual, False), (l0, None, False)]   # O_1m coincident : v1 site = v0 (shared index)


def op_vspec_point(x1, x2, letters, nsite):
    # POINT sigma^2: the two bilinear vertices PINNED to fixed sites x1,x2 (NO sum, NO A_x Y00).  Implemented
    # as one-hot (delta) weight vectors on two free site-letters -- perm_contrib_folded sums the letter, the
    # delta collapses the sum to the fixed site.  Machinery (make_config/AblkS/perm_contrib_folded) unchanged.
    # x1==x2 is the coincident-point case (carries the AblkS -1/2 contact, like O_1m).
    l0, l1 = letters
    e1 = np.zeros(nsite)
    e1[x1] = 1.0
    e2 = np.zeros(nsite)
    e2[x2] = 1.0
    return [(l0, e1, False), (l1, e2, False)]


def five_vertices():
    # the degree-5 (icosahedral) vertices: sites with 5 nearest neighbours (L1: all 12; L2: 12 of 42).
    # Uses the geom_hopping nns (same source as the stress-tensor agent's t00_ham_pole).
    import geom_hopping_claude as gh
    _om, _alpha, nns, nsite = gh.build(dc.GEOM, dc.L)
    return [i for i in range(nsite) if len(nns[i]) == 5]


def perm_contrib_folded(cycles, vtime, bA, vspec, Pmap):
    # kernels folded into A-blocks / weight vectors (no dense kernel matrix); s0 batch 'z' summed
    arrs = []
    idxs = []
    for cyc in cycles:
        m = len(cyc)
        for a in range(m):
            va = cyc[a]
            vb = cyc[(a + 1) % m]
            arr = bA[(vtime[va], vtime[vb])]
            if vspec[va][2]:
                arr = arr[:, Pmap, :, :, :]           # out-site (axis 1) -> antipode
            if vspec[vb][2]:
                arr = arr[:, :, :, Pmap, :]            # in-site (axis 3) -> antipode
            arrs.append(arr)
            idxs.append('z' + vspec[va][0] + SPIN[va] + vspec[vb][0] + SPIN[vb])
    for v in range(4):
        wv = vspec[v][1]
        if wv is not None:
            arrs.append(wv)
            idxs.append(vspec[v][0])
    sub = ','.join(idxs) + '->'
    path = _PATH_CACHE.get(sub)
    if path is None:
        path = np.einsum_path(sub, *arrs, optimize='optimal')[0]
        _PATH_CACHE[sub] = path
    val = np.einsum(sub, *arrs, optimize=path)
    return ((-1.0) ** len(cycles)) * val


def perm_contrib_batched(cycles, vtime, bA, Ksink, Ksrc):
    # bA[(c1,c2)] = batched A-block (nz, nsite,2,nsite,2) over the s0 batch ; sum over z (s0)
    arrs = []
    idxs = []
    for cyc in cycles:
        m = len(cyc)
        for a in range(m):
            va = cyc[a]
            vb = cyc[(a + 1) % m]
            arrs.append(bA[(vtime[va], vtime[vb])])
            idxs.append('z' + SITE[va] + SPIN[va] + SITE[vb] + SPIN[vb])
    arrs.append(Ksink)
    idxs.append('ij')
    arrs.append(Ksrc)
    idxs.append('kl')
    sub = ','.join(idxs) + '->'
    path = _PATH_CACHE.get(sub)
    if path is None:
        path = np.einsum_path(sub, *arrs, optimize='optimal')[0]
        _PATH_CACHE[sub] = path
    val = np.einsum(sub, *arrs, optimize=path)
    return ((-1.0) ** len(cycles)) * val


def matrix_one_config(k, KER, OFF, dual):
    # full 3x3 connected matrix C[i,j,dt] for FS (S-part + Stilde-part); s0 (source time) batched into einsum
    AblkS, AblkSt, twin, nsite, U, tau = make_config(k)
    d = SPLIT
    nop = len(KER)
    C = np.full((nop, nop, DTMAX), np.nan)
    omax = max(max(o) for o in OFF)
    dual = dual.astype(float)
    wY = dual * dc.Y00
    Pmap = antipodal_map()
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt + omax < twin and s + omax < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        offsets = set()
        for a in range(nop):
            for b in range(nop):
                vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                for va in range(4):
                    for vb in range(4):
                        offsets.add((vt[va], vt[vb]))
        bAS = {}
        bASt = {}
        for (c1, c2) in offsets:
            bAS[(c1, c2)] = np.array([AblkS(s + c1, s + c2) for s in s0s])
            if FS_STILDE_FURNISH_ENABLED:
                bASt[(c1, c2)] = np.array([AblkSt(s + c1, s + c2) for s in s0s])   # only the massive path calls AblkSt
        # FS S-part + Stilde-part.  DEFAULT (FS_STILDE_FURNISH_ENABLED=False, massless): the Stilde leg
        #   collapses to the S leg at m=0 (GW: S~ D_ov^{-dag} = tau, sigma_FS == sigma_PS), so AblkSt is NEVER
        #   CALLED (guard above) and the sum is 2 x the S-part.  The old code summed (bAS, bASt) with
        #   bASt = AblkSt = the buggy forward-furnished tau_gw leg (spurious bare D_ov^dag = an O(1) artifact).
        #   See fs_gw_collapse_v_agent_two_meson_claude.md.
        # >>> COMMAND: DO NOT DELETE THIS COMMENT OR THE FS_STILDE_FURNISH_ENABLED GUARD.  The guard disables
        #     the buggy tau_gw Stilde leg (AblkSt) at its call site; removing it re-introduces a fixed bug. <<<
        legs = (bAS, bASt) if FS_STILDE_FURNISH_ENABLED else (bAS, bAS)
        for a in range(nop):
            for b in range(nop):
                vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                vspec = op_vspec(a, ('i', 'j'), dual, wY) + op_vspec(b, ('k', 'l'), dual, wY)
                v = 0.0
                for bA in legs:             # Stilde==S at m=0 (GW collapse) -> 2 x S-part when disabled
                    for cyc in PERMS:
                        v += perm_contrib_folded(cyc, vt, bA, vspec, Pmap)
                C[a, b, dt] = v.real / len(s0s)
    return C


_WK = {}


def _init_worker(KER, OFF, dual):
    _WK["KER"] = KER
    _WK["OFF"] = OFF
    _WK["dual"] = dual


def _worker_matrix(k):
    return matrix_one_config(k, _WK["KER"], _WK["OFF"], _WK["dual"])


def gevp(Ct, C0):
    C0 = 0.5 * (C0 + C0.T)
    w, Uv = np.linalg.eigh(C0)
    keep = w > 1e-10 * w.max()
    Uk = Uv[:, keep] / np.sqrt(w[keep])
    M = Uk.T @ (0.5 * (Ct + Ct.T)) @ Uk
    return np.sort(np.linalg.eigvals(M).real)[::-1]


def run_gevp():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    nsite = dual.shape[0]
    Pmap = antipodal_map()
    d = SPLIT
    Kw = np.outer(dual * dc.Y00, dual * dc.Y00)
    K2m = np.zeros((nsite, nsite))
    for i in range(nsite):
        K2m[i, Pmap[i]] = dual[i]
    K1m = np.diag(dual.astype(float))
    KER = [Kw, K2m, K1m]
    OFF = [(0, 0), (0, 0), (0, d)]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    CACHEDIR = os.environ.get("CACHEDIR", "fs_gevp_cache_claude")   # in the working dir, not /tmp
    os.makedirs(CACHEDIR, exist_ok=True)

    if TIMING:
        import time
        t0 = time.time()
        C = matrix_one_config(ks[0], KER, OFF, dual)
        print("# TIMING one-config full 3x3: %.2f s  -> est %d cfg = %.1f min"
              % (time.time() - t0, len(dc.KS), (time.time() - t0) * len(dc.KS) / 60.0))
        for dt in range(1, 8):
            print("#  dt=%d  C[0,0]=% .3e C[1,1]=% .3e C[2,2]=% .3e C[0,1]=% .3e C[0,2]=% .3e C[1,2]=% .3e"
                  % (dt, C[0, 0, dt], C[1, 1, dt], C[2, 2, dt], C[0, 1, dt], C[0, 2, dt], C[1, 2, dt]))
        return

    tag = dc.ENS.split("nu0")[0]
    cache = "%s/fs_gevp_point_%s_%dcfg_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), len(ks), SPLIT)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded cached matrices <- %s  (%d cfg)" % (cache, allC.shape[0]))
    else:
        if NPROC > 1:
            import multiprocessing as mp
            print("# computing %d configs with %d workers ..." % (len(ks), NPROC))
            with mp.Pool(NPROC, initializer=_init_worker, initargs=(KER, OFF, dual)) as pool:
                res = pool.map(_worker_matrix, ks)
            allC = np.array(res)
        else:
            allC = np.array([matrix_one_config(k, KER, OFF, dual) for k in ks])   # (ncfg,3,3,DTMAX)
        np.save(cache, allC)
        print("# cached matrices -> %s" % cache)
    ncfg = allC.shape[0]
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))          # symmetrize
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    def effmass(Cm):
        ev = np.full((DTMAX, 3), np.nan)
        for dt in range(DTMAX):
            if np.any(~np.isfinite(Cm[:, :, dt])):
                continue
            try:
                ev[dt] = gevp(Cm[:, :, dt], Cm[:, :, T0])
            except Exception:
                pass
        with np.errstate(all="ignore"):
            return np.log(ev[:-1] / ev[1:])

    em_c = effmass(blk.mean(0))
    ems = np.array([effmass(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))

    print("\n#  t |   m0(err)        m1(err)        m2(err)")
    for t in range(T0, min(DTMAX - 1, 22)):
        print("#  %2d | %7.4f(%.4f)  %7.4f(%.4f)  %7.4f(%.4f)"
              % (t, em_c[t, 0], em_err[t, 0], em_c[t, 1], em_err[t, 1], em_c[t, 2], em_err[t, 2]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    M2PS = 0.644
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.7)
    ax.text(ts[-1] * 0.6, M2PS + 0.012, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    labs = ["state 0", "state 1", "state 2"]
    for n in range(3):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n], marker=mkr[n], ms=5, lw=1.1,
                    capsize=2.5, label=labs[n])
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(DTMAX - 1, 22))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FS connected GEVP {$\sigma^2_{00}$, $O_{2m}$(antipodal), $O_{1m}$(split)}  T0=%d  %s L1 %d cfg"
                 % (T0, tag, ncfg), fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_gevp_point_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    nsite = dual.shape[0]
    Pmap = antipodal_map()
    d = SPLIT
    # kernels (nsite x nsite) and time-offset pairs
    Kw = np.outer(dual * dc.Y00, dual * dc.Y00)                  # sigma^2_00
    K2m = np.zeros((nsite, nsite))
    for i in range(nsite):
        K2m[i, Pmap[i]] = dual[i]                                # O_2m antipodal
    K1m = np.diag(dual.astype(float))                            # O_1m coincident
    KER = [Kw, K2m, K1m]
    OFF = [(0, 0), (0, 0), (0, d)]                              # (sink/source) time offsets per op

    if VALIDATE:
        # C[0][0] (Y00) connected  vs  diag_effmass connected A,B,C,D,E,G sum
        de.CONTACT = 0.5
        k = dc.KS[0]
        AblkS, AblkSt, twin, nsite, U, tau = make_config(k)
        w00 = np.repeat(dual, NS) * dc.Y00
        Phi = [(U[a].conj().T) @ (w00[:, None] * U[a]) for a in range(twin)]
        CONN_IDX = [0, 1, 2, 3, 4, 6]
        print("# VALIDATION  C[0][0] (Y00) connected  vs  diags_pair CONN sum   (S-part only, leg=tau)")
        for dtc in (2, 4, 6):
            S = 0
            cnt = 0
            for s0 in range(twin):
                T = s0 + dtc
                if T >= twin:
                    continue
                tvec = [T, T, s0, s0]
                acc = 0.0
                for cyc in PERMS:
                    acc += perm_contrib(cyc, tvec, AblkS, Kw, Kw)
                S += acc
                cnt += 1
            S = (S / cnt).real
            dp = np.mean([(dc.W10 * de.diags_pair(Phi, tau, s0, s0 + dtc))[CONN_IDX].sum()
                          for s0 in range(twin - dtc)])
            print("#   dt=%d  point=% .8e  diags_CONN=% .8e  ratio=%.8f" % (dtc, S, dp, S / dp))
        return
    run_gevp()


if __name__ == "__main__":
    main()
