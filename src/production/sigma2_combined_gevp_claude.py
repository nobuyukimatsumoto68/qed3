#!/usr/bin/env python3
# sigma2_combined_gevp_claude.py -- ONE GEVP over the (2,2) SHELL single-meson operators AND all sigma^2 two-meson
#   operators, to resolve "the three states around 2m_PS" (the (2,2)~0.65, two-meson~0.72, excited~0.80) together.
#   Basis: {shell ell1/2,ell3/2,ell5/2 (PS single-meson, M-eigenspace windows)} + {sigma^2 PP: s2,O2m,O1m}.  NKEEP=4.
#   Blocks (per config, avg over source s, sink t=s+dt):
#     shell x shell : -Tr[P_shell(t) tau(t,s) P_shell(s) tau(s,t)]  (mode space).
#     sigma^2 x sigma^2 (PP) : flavfac-PP 4-vertex (perm_contrib_folded + FFPP), MODE_CONTACT.
#     CROSS shell x sigma^2 : 3-vertex triangle via the shell BLOB B(s,s)=A(s,t) K_shell(t) A(t,s),
#       K_shell(t)=U(t) P_shell(t) U(t)^H ; cross_op = -Tr[B Wa A(s,s)_sub Wb] (+a<->b) with op geometry (Wa/Wb/antipode).
#   VALIDATE=1: with K_shell = diag(wY) (local psibar psi sink) the s2 cross must match sigma2_mPS_gevp C_12.
#   Run: MODE_CONTACT=1 ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 WINDOWS=0-4,4-12,12-24 OPS2=0,1,2 NKEEP=4 REBT=3 T0=2 \
#        BINSIZE=10 NCFG=200 python3 sigma2_combined_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MODE_CONTACT", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import sigma2_mPS_gevp_claude as MG
import hankel_rebase_scan_claude as hs

WINDOWS = [tuple(int(x) for x in w.split("-")) for w in os.environ.get("WINDOWS", "0-4,4-12,12-24").split(",")]
OPS2 = [int(x) for x in os.environ.get("OPS2", "0,1,2").split(",")]   # sigma^2 geometries (0=s2,1=O2m,2=O1m)
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
REBT = int(os.environ.get("REBT", "3"))
T0 = int(os.environ.get("T0", "2"))
NKEEP = int(os.environ.get("NKEEP", "4"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0").split(",")]
VALIDATE = int(os.environ.get("VALIDATE", "0"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))
NS = G.NS


def shell_projector(M, lo, hi):
    w, U = np.linalg.eigh(1j * M)
    sel = np.argsort(-np.abs(w))[lo:hi]
    Us = U[:, sel]
    return Us @ Us.conj().T


def antipode_perm(nsite, Pmap):
    # site-antipode permutation lifted to the 2Ns (site,spin) flat index j = NS*x + s
    p = np.empty(nsite * NS, dtype=int)
    for x in range(nsite):
        for s in range(NS):
            p[NS * x + s] = NS * Pmap[x] + s
    return p


def op_src_middle(op, A00_flat, wY_flat, dual_flat, Pap):
    # the sigma^2 SOURCE "middle" Wa A(s,s)_sub Wb as a (2Ns,2Ns) matrix, per geometry (op_vspec structure).
    # op0 s2: diag(wY) A00 diag(wY) ; op1 O2m: diag(dual) A00[antipode cols] (2nd vertex None-weight at P(site)) ;
    # op2 O1m: diag(dual) A00 (2nd vertex None-weight, coincident).
    if op == 0:
        return (wY_flat[:, None] * A00_flat) * wY_flat[None, :]
    if op == 1:
        return (dual_flat[:, None] * A00_flat)[:, Pap]           # 2nd vertex at antipode (weight 1), col-reindex
    return (dual_flat[:, None] * A00_flat)                        # O1m coincident (weight 1)


def one_config(k, dualf, wY, Pmap):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    twoNs = nsite * NS
    Pap = antipode_perm(nsite, Pmap)
    wY_flat = np.repeat(dualf, NS) * dc.Y00           # (2Ns,) spinor vertex weight for s2 (area*Y00)
    dual_flat = np.repeat(dualf, NS)                  # (2Ns,) area weight for O2m/O1m first vertex
    nsh = len(WINDOWS)
    nop = len(OPS2)
    ntot = nsh + nop
    # shell projectors + nonlocal K_shell(t) = U P_shell U^H (flat 2Ns x 2Ns) per timeslice
    Ps = [[shell_projector(tau[a, a] - 0.5 * np.eye(tau.shape[-1]), lo, hi) for (lo, hi) in WINDOWS] for a in range(twin)]
    Ksh = [[U[a] @ Ps[a][w] @ U[a].conj().T for w in range(nsh)] for a in range(twin)]  # (2Ns,2Ns)
    C = np.full((ntot, ntot, DTMAX), np.nan)

    # sigma^2 x sigma^2 (PP): reuse the perm_contrib_folded machinery via sigma2_mPS one_config's structure
    vsp = {op: G.op_vspec(op, ('i', 'j'), dualf, wY) for op in OPS2}
    vsp_src = {op: G.op_vspec(op, ('k', 'l'), dualf, wY) for op in OPS2}

    for dt in range(DTMAX):
        s0s = np.array([s for s in range(twin) if s + dt < twin])
        if len(s0s) == 0:
            continue
        offs = {(0, 0), (dt, dt), (0, dt), (dt, 0)}
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        # --- shell x shell (mode space) ---
        for a in range(nsh):
            for b in range(nsh):
                acc = 0.0
                for s in s0s:
                    t = s + dt
                    acc += -np.trace(Ps[t][a] @ tau[t, s] @ Ps[s][b] @ tau[s, t]).real
                C[a, b, dt] = acc / len(s0s)
        # --- sigma^2 x sigma^2 (PP) 4-vertex ---
        vt = [dt, dt, 0, 0]
        for ia, opa in enumerate(OPS2):
            for ib, opb in enumerate(OPS2):
                vspec = G.op_vspec(opa, ('i', 'j'), dualf, wY) + G.op_vspec(opb, ('k', 'l'), dualf, wY)
                base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap) for cyc in G.PERMS])
                C[nsh + ia, nsh + ib, dt] = (MG.FFPP @ base).real / len(s0s)
        # --- CROSS shell x sigma^2 via blob B(s,s) = A(s,t) K_shell(t) A(t,s) ---
        for a in range(nsh):
            for ib, opb in enumerate(OPS2):
                acc = 0.0
                for s in s0s:
                    t = s + dt
                    Ats = (U[t] @ tau[t, s] @ U[s].conj().T)                # A(t,s) flat 2Ns
                    Ast = (U[s] @ tau[s, t] @ U[t].conj().T)                # A(s,t) flat 2Ns
                    A00 = (U[s] @ tau[s, s] @ U[s].conj().T) - 0.5 * np.eye(twoNs)   # A(s,s)_sub (MODE-space contact below)
                    if int(os.environ.get("MODE_CONTACT", "1")):
                        A00 = U[s] @ (tau[s, s] - 0.5 * np.eye(tau.shape[-1])) @ U[s].conj().T
                    B = Ast @ Ksh[t][a] @ Ats                               # blob at time s (2Ns,2Ns)
                    Msrc = op_src_middle(opb, A00, wY_flat, dual_flat, Pap)
                    acc += -2.0 * np.trace(B @ Msrc).real                   # (+a<->b factor 2 for s2; approx for O2m/O1m)
                cval = acc / len(s0s)
                C[a, nsh + ib, dt] = cval
                C[nsh + ib, a, dt] = cval
    return C


def gevp(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    Vp = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, Vp, T0), Vp


def main():
    tag = dc.ENS.split("nu0")[0]
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d  COMBINED shell+sigma^2 GEVP  windows=%s ops2=%s  %d cfg  reb%d@%d T0=%d NKEEP=%d MC=%s"
          % (tag, dc.L, WINDOWS, OPS2, len(ks), NKEEP, REBT, T0, NKEEP, os.environ.get("MODE_CONTACT")))
    allC = np.array([one_config(k, dualf, wY, Pmap) for k in ks])
    ncfg = allC.shape[0]
    nb = max(ncfg // BINSIZE, 1)
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = gevp(blk.mean(0), None)
    if nb < 2:
        em_e = np.zeros_like(em_c)
    else:
        ems = np.array([gevp(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
        em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax, nk = em_c.shape
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(nk)))
    for t in range(T0, min(tmax, 16)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(nk))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.7)
    ax.text(TMAXPLOT * 0.7, m2ps + 0.01, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="firebrick", ls=":", lw=1.2, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, mps + 0.01, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="firebrick")
    cols = ["firebrick", "tab:purple", "tab:blue", "tab:green", "tab:orange"]
    mkr = ["o", "D", "s", "^", "v"]
    for n in range(nk):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 5], marker=mkr[n % 5], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.set_ylim(0.2, 1.05)
    ax.set_xlim(T0, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"COMBINED (2,2)-shell + $\sigma^2$ GEVP  %s L%d %dcfg  win=%s ops2=%s NKEEP=%d"
                 % (tag, dc.L, ncfg, WINDOWS, OPS2, NKEEP), fontsize=9)
    ax.legend(fontsize=9, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_combined_gevp_%s_L%d_claude.png" % (tag, dc.L)
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
