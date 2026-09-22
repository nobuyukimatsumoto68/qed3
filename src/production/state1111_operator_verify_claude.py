#!/usr/bin/env python3
# state1111_operator_verify_claude.py
#   Verify the {1,1,1,1}=(2,2) operator picture in the FREE limit.
#
#   NM's picture: the (2,2) one-meson operator is an EQUAL-TIME nonlocal bilinear
#       O(t) = P_{ell=0} [ psibar(t,x) Xi_{3/2}(x,y) Xi_{3/2}(y,z) psi(t,z) ]
#   where Xi_{ell} := D^{-1}|_{ell} is the SPATIAL propagator restricted to the ell-shell (Xi = the 2D
#   spinor basis, qed3_v2-6.pdf Eq C.18), NOT a bare projector.  BOTH fermion legs are dressed to ell=3/2,
#   so the state is a single meson with both constituents in the ell=3/2 shell -> energy 2 E_{3/2}.
#   On a single (degenerate) shell, Xi_ell^2 = (1/lambda^2)(P_+ + P_-) is proportional to the shell density
#   P_ell = sum_{m,i3} Xi_a Xi_a^dag, so the ell=0 scalar operator O_ell = psibar P_ell psi restricts BOTH
#   legs to the shell.  Two-point with a spatial kernel K reduces to the standard distillation trace
#       <O(t) O(0)> = -Tr[ Phi_K(t) tau(t,s) Phi_K(s) tau(s,t) ] ,   Phi_K = V^dag K V .
#
#   COMPLICATION (found): the continuum C.18 shells are exactly orthonormal on the grid, but the free
#   *lattice* D_ov has only ICOSAHEDRAL symmetry (not full SO(3)), so tau MIXES continuum shells -> a single
#   P_ell operator is ground-contaminated.  A GEVP of shell operators {P_{1/2}, P_{3/2}, ...} separates them:
#   state0 -> m_PS (both legs ell=1/2), state1 -> (2,2) (both legs ell=3/2).  The mixing shrinks with
#   refinement (L1 12 sites obstructs state1; L2 42 sites recovers it near 2E_{3/2}).
#
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 state1111_operator_verify_claude.py
#         ENS=free LREF=2 NVDIR=distill_Nv84 python3 state1111_operator_verify_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "1")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import scipy.linalg as sla
import distill_contract_claude as dc
import free_wavefunctions_claude as fw

T0 = int(os.environ.get("T0", "2"))
# references (free): m_PS = 2 E_{1/2}, (2,2) = 2 E_{3/2}
M_PS = {1: 0.378, 2: 0.393}[dc.L]
STATE22 = {1: 0.556, 2: 0.690}[dc.L]
SHELLS = {1: [1, 2], 2: [1, 2, 3]}[dc.L]     # lambda-shells in the GEVP basis


def rot_and_su2(axis, angle):
    # 3x3 rotation R (Rodrigues) and its 2x2 SU(2) rep U (spinor rotation)
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    U = np.cos(angle / 2) * np.eye(2) - 1j * np.sin(angle / 2) * (n[0] * sx + n[1] * sy + n[2] * sz)
    return R, U


def build_shell_Xi(lam, sites, R, U):
    # shell modes with a ROTATED quantization axis so no lattice site hits a pole:
    #   Xi^{(R)}(x) = U . Xi^std(R^{-1} x)  (U = SU(2) spinor rotation).  P_ell = sum |Xi><Xi| is
    #   rotation-invariant, so this is the projector at the fixed sites x.
    xr = sites @ R                                     # R^{-1} x = R^T x  (rows)
    zr = np.clip(xr[:, 2], -1.0, 1.0)
    thr = np.arccos(zr)
    phr = np.mod(np.arctan2(xr[:, 1], xr[:, 0]), 2.0 * np.pi)
    assert np.all(np.abs(np.abs(zr) - 1.0) > 1e-6), "a site still at a pole after rotation"
    modes = fw.shell_modes(lam)
    cols = []
    for (m, n, i3) in modes:
        up, dn = fw.psi(m, n, i3, thr, phr)
        sp = np.stack([up, dn], axis=1) @ U.T          # apply spinor rotation per site
        cols.append(sp.ravel())                        # (2Ns,) site-major spin-minor
    return np.stack(cols, axis=1), modes


def shell_vertex(V_t, Xi, area2):
    # Phi_ell(t) = V^dag P_ell V = O O^H,  O[c,a] = <V_c | Xi_a> (area-weighted)
    O = (area2[None, :] * V_t.conj()) @ Xi
    return O @ O.conj().T


def main():
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    dual = dc.dual_areas_from_mesh()
    area2 = np.repeat(dual, dc.NS)
    Rrot, Urot = rot_and_su2([1.0, 2.0, 3.0], 0.9)
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = tau.shape[-1]
    print("# FREE L%d  {1,1,1,1}=(2,2) operator verify (shell GEVP):  Nv=%d nsite=%d twin=%d" % (dc.L, Nv, len(dual), twin))
    print("# operators = psibar P_ell psi (both legs in shell ell), P_ell from Eq C.18 modes ; GEVP T0=%d" % T0)
    print("# PREDICT: state0 = m_PS = %.3f (ell=1/2) ; state1 = (2,2) = 2E_{3/2} = %.3f (ell=3/2)\n" % (M_PS, STATE22))

    Phi = {}
    for lam in SHELLS:
        Xi, modes = build_shell_Xi(lam, sites, Rrot, Urot)
        G = (area2[:, None] * Xi).conj().T @ Xi
        offd = np.linalg.norm(G - np.diag(np.diag(G))) / np.linalg.norm(np.diag(G))
        Phi[lam] = [shell_vertex(V[tsrc0 + a].T, Xi, area2) for a in range(twin)]
        print("#   shell lambda=%d (ell=%.1f)  n=%d modes  mode-overlap ||offdiag||/||diag||=%.1e"
              % (lam, lam - 0.5, len(modes), offd))

    nop = len(SHELLS)
    Cmat = np.zeros((nop, nop, twin))
    for dt in range(twin):
        ns = twin - dt
        for i, li in enumerate(SHELLS):
            for j, lj in enumerate(SHELLS):
                acc = 0.0
                for s in range(ns):
                    t = s + dt
                    acc += (-np.trace(Phi[li][t] @ tau[t, s] @ Phi[lj][s] @ tau[s, t])).real
                Cmat[i, j, dt] = acc / ns

    C0 = 0.5 * (Cmat[:, :, T0] + Cmat[:, :, T0].T)
    lam_t = np.full((twin, nop), np.nan)
    for dt in range(twin):
        Ct = 0.5 * (Cmat[:, :, dt] + Cmat[:, :, dt].T)
        try:
            lam_t[dt] = np.sort(sla.eigvals(Ct, C0).real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam_t[:-1] / lam_t[1:])

    print("\n#  t | " + "  ".join("state%d" % n for n in range(nop)) + "     (m_PS=%.3f  (2,2)=%.3f)" % (M_PS, STATE22))
    for dt in range(1, min(20, twin - 1)):
        print("#  %2d |  %s" % (dt, "  ".join("%7.4f" % em[dt, n] if np.isfinite(em[dt, n]) else "   --- " for n in range(nop))))

    # plateau windows: state0 late (ground); state1 mid (before contamination)
    def med(a, lo, hi):
        w = a[lo:hi]
        w = w[np.isfinite(w)]
        return np.median(w) if w.size else np.nan
    e0 = med(em[:, 0], 8, 14)
    e1 = med(em[:, 1], 6, 10)
    print("\n# state0 (ell=1/2) plateau[8-14] = %.4f   [m_PS=%.3f]   %s" % (e0, M_PS, "OK" if abs(e0 - M_PS) < 0.02 else "check"))
    print("# state1 (ell=3/2) plateau[6-10] = %.4f   [(2,2)=%.3f]   %s" % (e1, STATE22, "OK/near" if abs(e1 - STATE22) < 0.06 else "obstructed (shell-mix)"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(twin - 1)
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    ax.axhline(M_PS, color="tab:green", ls="--", lw=1, alpha=0.6)
    ax.text(twin * 0.55, M_PS + 0.006, r"$m_{PS}=2E_{1/2}=%.3f$" % M_PS, color="tab:green", fontsize=9)
    ax.axhline(STATE22, color="tab:red", ls="--", lw=1, alpha=0.6)
    ax.text(twin * 0.55, STATE22 + 0.006, r"$(2,2)=2E_{3/2}=%.3f$" % STATE22, color="tab:red", fontsize=9)
    cols = ["tab:green", "tab:red", "tab:gray"]
    mk = ["s", "o", "^"]
    for n in range(nop):
        g = np.isfinite(em[:, n])
        ax.plot(ts[g], em[g, n], color=cols[n % 3], marker=mk[n % 3], ms=5, lw=1.1, label="state %d" % n)
    ax.set_ylim(0.25, max(1.0, STATE22 + 0.2))
    ax.set_xlim(0, min(twin - 2, 18))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L%d  shell-operator GEVP: $\bar\psi P_\ell\psi$ (both legs in shell $\ell$)" % dc.L)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/state1111_operator_verify_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
