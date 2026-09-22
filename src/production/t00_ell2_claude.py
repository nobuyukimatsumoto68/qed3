#!/usr/bin/env python3
# t00_ell2_claude.py
# Stress tensor in the CORRECT channel: ell=2 projection of the e.sigma energy density O_H (plan: t00_ell2_impl_plan_claude.md).
#   O_H^(m)(t) = eta^H W^(m)(t) xi + h.c.,   W^(m)_ij = W_ij * Y_lm(n_mid(ij)),  REAL Y_lm  (W^(m) stays anti-hermitian)
#   => Phi_WH = -Phi_W, the two halves of the connected 2pt are equal:
#   C_conn(dt) = -2 Re < Tr[Phi_m(t) tau(t,s) Phi_m(s) tau(s,t)] >_{s,m}        (dt>=1; dt=0 contact not treated)
#   connected-only = flavor-ADJOINT spin-2 (Delta = 3+gamma_adj); singlet needs the disconnected loops:
#   L_m(t) = -2 Re Tr[Phi_m(t) (tau(t,t) - 1/2)]   (W half leg tau, W^dag half leg 1-tau)  -> cached per config.
# ell=0 is kept as a reference channel ((n,n) pair states only: 2E_{1/2}, 2E_{3/2}, ...).
# Run free L1:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_ell2_claude.py
# Run free L2:  ENS=free NVDIR=distill_Nv84 LREF=2 python3 t00_ell2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import geom_hopping_claude as gh
import t00_stress_ham_interacting_claude as ti

NS = dc.NS
LREF = int(os.environ.get("LREF", "1"))
GEOM = os.environ.get("GEOM", "../../geometry/data/")
ELLS = [int(x) for x in os.environ.get("ELLS", "0,2").split(",")]
DTMAX = int(os.environ.get("DTMAX", "24"))
KMIN = int(os.environ.get("KMIN", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
FREE = (dc.ENS == "free")
CONFIGDIR = os.environ.get("CONFIGDIR", dc.ENS)
ROVER = 1.0 / 0.189


def real_ylm(n, ell):
    # real spherical harmonics (Cartesian forms) at unit vectors n (npts,3) -> (2ell+1, npts)
    x = n[:, 0]
    y = n[:, 1]
    z = n[:, 2]
    if ell == 0:
        return np.array([np.full(len(x), 0.5 / np.sqrt(np.pi))])
    if ell == 1:
        c = np.sqrt(3.0 / (4.0 * np.pi))
        return np.array([c * x, c * y, c * z])
    if ell == 2:
        c = 0.5 * np.sqrt(15.0 / np.pi)
        return np.array([c * x * y,
                         c * y * z,
                         c * z * x,
                         0.5 * c * (x * x - y * y),
                         0.25 * np.sqrt(5.0 / np.pi) * (3.0 * z * z - 1.0)])
    raise ValueError("ell must be 0,1,2")


def link_weight_matrices(pts, nns, nsite, ell):
    # Ymat[m] (nsite,nsite): Y_lm at the normalized midpoint of link (i,j); symmetric in (i,j)
    nm = 2 * ell + 1
    Ymat = np.zeros((nm, nsite, nsite))
    for i in range(nsite):
        for j in nns[i]:
            mid = pts[i] + pts[j]
            mid = mid / np.linalg.norm(mid)
            Ymat[:, i, j] = real_ylm(mid[None, :], ell)[:, 0]
    return Ymat


def corr_one(k, geo, Ymats):
    om, alpha, nns, nsite, tab = geo
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    n_links = len(tab) // 2
    if FREE:
        sp = np.zeros((V.shape[0], n_links))
    else:
        sp = ti.read_config_sp("%s/ckpoint_lat.%d" % (CONFIGDIR, k), n_links)
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    dtmax = min(DTMAX, twin)
    out_c = {}
    out_l = {}
    Wt = []
    Vts = []
    for a in range(twin):
        t = tsrc0 + a
        Vts.append(V[t].T)
        Wt.append(ti.build_W_gauge(om, alpha, nns, nsite, tab, sp[t]))     # r=0 anti-hermitian hop
    for ell in ELLS:
        Ymat = Ymats[ell]
        nm = Ymat.shape[0]
        Phi = np.zeros((twin, nm, Nv, Nv), complex)
        for a in range(twin):
            for m in range(nm):
                Wm = Wt[a] * np.kron(Ymat[m], np.ones((NS, NS)))
                Phi[a, m] = Vts[a].conj().T @ Wm @ Vts[a]
        C = np.zeros(dtmax)
        for dt in range(1, dtmax):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                M = Phi[t] @ tau[t, s]
                N = Phi[s] @ tau[s, t]
                acc += np.einsum("mpq,mqp->", M, N).real
            C[dt] = -2.0 * acc / (ns * nm)
        L = np.zeros((twin, nm))
        for a in range(twin):
            G = tau[a, a] - 0.5 * Iv
            L[a] = -2.0 * np.einsum("mpq,qp->m", Phi[a], G).real
        out_c[ell] = C
        out_l[ell] = L
    return out_c, out_l


def main():
    om, alpha, nns, nsite = gh.build(GEOM, LREF)
    tab = ti.load_link_table("primal_links_n%d_claude.dat" % LREF)
    pts = dc.load_vec3(GEOM + "pts_n%d.dat" % LREF)
    geo = (om, alpha, nns, nsite, tab)
    Ymats = {ell: link_weight_matrices(pts, nns, nsite, ell) for ell in ELLS}
    ks = dc.KS if FREE else [k for k in dc.KS if k >= KMIN]
    if NCFG:
        ks = ks[:NCFG]
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s NVDIR=%s LREF=%d nsite=%d ncfg=%d ELLS=%s" % (tag, dc.NVDIR, LREF, nsite, len(ks), ELLS))
    cdir = os.environ.get("CACHEDIR", "t00_ham_cache_claude")
    os.makedirs(cdir, exist_ok=True)
    cache = "%s/ell2_%s_%s_n%d_km%d_claude.npz" % (cdir, tag.replace(".", "p"), dc.NVDIR, len(ks), KMIN)
    if os.path.exists(cache):
        d = np.load(cache)
        allC = {ell: d["C%d" % ell] for ell in ELLS}
        print("# loaded cache <- %s" % cache)
    else:
        allC = {ell: [] for ell in ELLS}
        allL = {ell: [] for ell in ELLS}
        for i, k in enumerate(ks):
            c, l = corr_one(k, geo, Ymats)
            for ell in ELLS:
                allC[ell].append(c[ell])
                allL[ell].append(l[ell])
            if (i + 1) % 50 == 0:
                print("#   ... %d/%d" % (i + 1, len(ks)))
        save = {}
        for ell in ELLS:
            allC[ell] = np.array(allC[ell])
            save["C%d" % ell] = allC[ell]
            save["L%d" % ell] = np.array(allL[ell])
        save["ks"] = np.array(ks)
        np.savez(cache, **save)
        print("# cached -> %s" % cache)

    em = {}
    er = {}
    for ell in ELLS:
        A = allC[ell]
        ncfg = A.shape[0]
        nb = max(ncfg // BINSIZE, 1)
        if nb > 1:
            blk = np.array([A[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
        else:
            blk = A.mean(0)[None, :]
        with np.errstate(all="ignore"):
            em[ell] = np.log(blk.mean(0)[:-1] / blk.mean(0)[1:])
            if nb > 1:
                jm = np.array([np.log(np.delete(blk, i, 0).mean(0)[:-1] / np.delete(blk, i, 0).mean(0)[1:]) for i in range(nb)])
                er[ell] = np.sqrt((nb - 1) * np.nanmean((jm - np.nanmean(jm, 0)) ** 2, axis=0))
            else:
                er[ell] = np.zeros_like(em[ell])

    print("\n# refs: 2/R=%.3f  3/R=%.3f  4/R=%.3f" % (2 / ROVER, 3 / ROVER, 4 / ROVER))
    hdr = "#  dt |"
    for ell in ELLS:
        hdr += "  ell=%d m_eff(err)      C(dt)     |" % ell
    print(hdr)
    for dt in range(1, DTMAX - 1):
        row = "#  %2d |" % dt
        for ell in ELLS:
            if dt < em[ell].shape[0]:
                row += "  %8.4f(%.4f)  %11.3e |" % (em[ell][dt], er[ell][dt], allC[ell].mean(0)[dt])
        print(row)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sty = {0: ("tab:red", "o"), 1: ("tab:green", "^"), 2: ("tab:blue", "s")}
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for y, lab in [(2 / ROVER, "2/R"), (3 / ROVER, "3/R"), (4 / ROVER, "4/R")]:
        ax.axhline(y, color="gray", ls="--", lw=0.9)
        ax.text(DTMAX - 3.5, y + 0.006, lab, fontsize=9, color="gray")
    for ell in ELLS:
        ts = np.arange(em[ell].shape[0])
        g = np.isfinite(em[ell]) & (ts >= 1) & (er[ell] < 0.3)
        col, mk = sty[ell]
        ax.errorbar(ts[g], em[ell][g], yerr=er[ell][g], color=col, marker=mk, ms=5, lw=1.0, capsize=2.5,
                    label=r"$O_H$ $\ell=%d$ (connected)" % ell)
    ax.set_ylim(0.2, 1.2)
    ax.set_xlim(0, DTMAX - 1)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$O_H$ $\ell$-projected  %s  %s  L%d" % (tag, dc.NVDIR, LREF), fontsize=10.5)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_ell2_%s_L%d_claude.png" % (tag, LREF)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
