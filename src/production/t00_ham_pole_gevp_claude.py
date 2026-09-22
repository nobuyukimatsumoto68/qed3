#!/usr/bin/env python3
# t00_ham_pole_gevp_claude.py
# POINT (local) energy-density two-point: put the e.sigma Hamiltonian density at the NORTH and SOUTH pole sites
# (the sum of hoppings in each pole's star), NOT the Y00/global sum, and build the 2x2 correlator matrix + GEVP.
# On the sphere the symmetric N+S combination is ell=0 (energy density), antisymmetric N-S is ell=1.
# Vertex: W_pole(x0) = the star of site x0 (blocks of the full gauge-dressed Wilson W touching x0), r=0.
# Contraction (elemental; same loop as the global O_H): C_ij(dt) = -(1/N) sum_a { Tr[Phi_i(b) tau(b,a) Phi_j(a) tau(a,b)]
#   + Tr[Phi_iH(b) tau(b,a) Phi_jH(a) tau(a,b)] },  Phi_i(a)=V(a)^H W_i(a) V(a),  i,j in {N,S},  b=a+dt.
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000 NVDIR=distill_Nv24_v2 \
#       LREF=2 python3 t00_ham_pole_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as fg
import geom_hopping_claude as gh
import t00_stress_ham_interacting_claude as ti

NS = dc.NS
DTMAX = int(os.environ.get("DTMAX", "24"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
KMIN = int(os.environ.get("KMIN", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
LREF = int(os.environ.get("LREF", "1"))
GEOM = os.environ.get("GEOM", "../../geometry/data/")
CONFIGDIR = os.environ.get("CONFIGDIR", dc.ENS)
T0 = int(os.environ.get("T0", "3"))
ROVER = 1.0 / 0.189


def pole_sites():
    pts = dc.load_vec3(GEOM + "pts_n%d.dat" % LREF)
    north = int(np.argmax(pts[:, 2]))
    south = int(np.argmin(pts[:, 2]))
    return north, south


def mask_star(W, x0, nsite):
    # keep only blocks of W touching site x0 (row x0 and col x0) = the star (local density at x0)
    Wl = np.zeros_like(W)
    Wl[NS * x0:NS * x0 + 2, :] = W[NS * x0:NS * x0 + 2, :]
    Wl[:, NS * x0:NS * x0 + 2] = W[:, NS * x0:NS * x0 + 2]
    return Wl


def five_vertices(nns, nsite):
    return [i for i in range(nsite) if len(nns[i]) == 5]     # the 12 icosahedral (degree-5) sites


def corr_matrix(k, om, alpha, nns, nsite, tab, opsites, add_global):
    # basis: a point-operator (star) at each site in opsites, plus optionally the global H (ell=0 full sum)
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    sp = ti.read_config_sp("%s/ckpoint_lat.%d" % (CONFIGDIR, k), len(tab) // 2)
    N = NS * nsite
    nop = len(opsites) + (1 if add_global else 0)
    Phi = np.zeros((twin, nop, V.shape[1], V.shape[1]), complex)     # (twin, nop, Nv, Nv)
    PhiH = np.zeros_like(Phi)
    for a in range(twin):
        Vt = V[tsrc0 + a].T
        W = ti.build_W_gauge(om, alpha, nns, nsite, tab, sp[tsrc0 + a])
        verts = [mask_star(W, x0, nsite) for x0 in opsites]
        if add_global:
            verts.append(W)
        for o in range(nop):
            Phi[a, o] = Vt.conj().T @ verts[o] @ Vt
            PhiH[a, o] = Vt.conj().T @ verts[o].conj().T @ Vt
    C = np.full((nop, nop, DTMAX), np.nan)
    for dt in range(DTMAX):
        na = twin - dt
        if na <= 0:
            continue
        acc = np.zeros((nop, nop), complex)
        for a in range(na):
            b = a + dt
            M = Phi[b] @ tau[b, a]           # (nop, Nv, Nv)
            Ncol = Phi[a] @ tau[a, b]
            MH = PhiH[b] @ tau[b, a]
            NH = PhiH[a] @ tau[a, b]
            acc += np.einsum("ipq,jqp->ij", M, Ncol) + np.einsum("ipq,jqp->ij", MH, NH)
        C[:, :, dt] = (-acc / na).real
    return C


def main():
    om, alpha, nns, nsite = gh.build(GEOM, LREF)
    tab = ti.load_link_table("primal_links_n%d_claude.dat" % LREF)
    opsites = five_vertices(nns, nsite)                       # 12 icosahedral degree-5 vertices
    add_global = int(os.environ.get("ADDGLOBAL", "1"))
    nop_tot = len(opsites) + (1 if add_global else 0)
    print("# ENS=%s NVDIR=%s nsite=%d  #5-vertices=%d  addH=%d -> nop=%d  CONFIGDIR=%s"
          % (dc.ENS, dc.NVDIR, nsite, len(opsites), add_global, nop_tot, CONFIGDIR))
    ks = [k for k in dc.KS if k >= KMIN]
    if NCFG:
        ks = ks[:NCFG]
    tag = dc.ENS.split("nu0")[0]
    cdir = os.environ.get("CACHEDIR", "t00_ham_cache_claude")
    os.makedirs(cdir, exist_ok=True)
    cache = "%s/pole12_%s_%s_n%d_dt%d_km%d_g%d_claude.npy" % (cdir, tag.replace(".", "p"), dc.NVDIR, len(ks), DTMAX, KMIN, add_global)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded cache <- %s (%d cfg)" % (cache, allC.shape[0]))
    else:
        allC = []
        for i, k in enumerate(ks):
            allC.append(corr_matrix(k, om, alpha, nns, nsite, tab, opsites, add_global))
            if (i + 1) % 50 == 0:
                print("#   ... %d/%d" % (i + 1, len(ks)))
        allC = np.array(allC)
        np.save(cache, allC)
        print("# cached -> %s" % cache)
    ncfg = allC.shape[0]
    nop = allC.shape[1]
    hidx = nop - 1 if add_global else 0                   # reference diagonal (global H if present, else a pole)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))          # symmetrize
    sgn = np.sign(np.nanmean(allC[:, hidx, hidx, 3]))
    allC = sgn * allC

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    def gevp_em(Cm):
        ev = np.full((DTMAX, nop), np.nan)
        for dt in range(DTMAX):
            if np.any(~np.isfinite(Cm[:, :, dt])) or np.any(~np.isfinite(Cm[:, :, T0])):
                continue
            try:
                e = fg.gevp(Cm[:, :, dt], Cm[:, :, T0])
                ev[dt, :len(e)] = e[:nop]
            except Exception:
                pass
        with np.errstate(all="ignore"):
            return np.log(ev[:-1] / ev[1:])

    em_c = gevp_em(blk.mean(0))
    ems = np.array([gevp_em(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, axis=0))

    # global-H diagonal point-to-point (single-op reference)
    Cm = blk.mean(0)
    with np.errstate(all="ignore"):
        emHH = np.log(Cm[hidx, hidx, :-1] / Cm[hidx, hidx, 1:])

    print("# ncfg=%d bin%d nb=%d T0=%d  (%dx%d GEVP: 12 five-vertices%s)"
          % (ncfg, BINSIZE, nb, T0, nop, nop, "+globalH" if add_global else ""))
    print("\n#  t | GEVP m0(err)    m1(err)    m2(err)   | ref point2point")
    for t in range(T0, min(em_c.shape[0], 20)):
        print("#  %2d | %7.4f(%.4f) %7.4f(%.4f) %7.4f(%.4f) | %7.4f"
              % (t, em_c[t, 0], em_err[t, 0], em_c[t, 1], em_err[t, 1], em_c[t, 2], em_err[t, 2],
                 emHH[t] if t < emHH.shape[0] else np.nan))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    AX = float(os.environ.get("AXIAL", "0"))
    ts = np.arange(em_c.shape[0])
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    ax.axhline(3.0 / ROVER, color="tab:red", ls="--", lw=1, alpha=0.6)
    ax.text(T0 + 0.1, 3.0 / ROVER + 0.008, r"free $T_{00}=0.567$", fontsize=9, color="tab:red")
    if AX:
        ax.axhline(AX, color="tab:purple", ls=":", lw=1.2, alpha=0.8)
        ax.text(T0 + 0.1, AX + 0.008, r"axial=%.3f" % AX, fontsize=9, color="tab:purple")
        ax.axhline(1.5 * AX, color="tab:green", ls="-.", lw=1.2, alpha=0.8)
        ax.text(T0 + 0.1, 1.5 * AX + 0.008, r"1.5xaxial=%.3f" % (1.5 * AX), fontsize=9, color="tab:green")
    g0 = np.isfinite(em_c[:, 0]) & (em_err[:, 0] < 0.4)
    ax.errorbar(ts[g0], em_c[g0, 0], yerr=em_err[g0, 0], color="tab:blue", marker="o", ms=5, lw=1.2, capsize=2.5,
                label="GEVP ground (ell=0)")
    g1 = np.isfinite(em_c[:, 1]) & (em_err[:, 1] < 0.4)
    ax.errorbar(ts[g1], em_c[g1, 1], yerr=em_err[g1, 1], color="tab:orange", marker="s", ms=4.5, lw=1.0, capsize=2.5,
                markerfacecolor="none", label="GEVP excited (ell=1)")
    gp = np.isfinite(emHH)
    ax.plot(np.arange(emHH.shape[0])[gp], emHH[gp], "-", color="lightgray", lw=1.0, marker=".", ms=4, label="global-H point2point")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlim(T0, min(em_c.shape[0], 20))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$T_{00}$ 12 five-vertices + global H GEVP  %s  %d cfg" % (tag, ncfg), fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_ham_pole_gevp_%s_claude.png" % tag
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
