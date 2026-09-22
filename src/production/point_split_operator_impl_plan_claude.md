# Point-split scalar operator for the $\ell=0$ two-meson GEVP

## Goal / physics

In the free $\ell=0$ $\sigma_{00}^2$ correlator, diagram A tracks the **first radially-excited $\ell=0$
single meson** $(2,2)$ = both fermions in the $\lambda=2$ Dirac shell (Eq C.17 of `qed3_v2-6.pdf`,
$\lambda_{|m|,n}=n+|m|+\tfrac12$), $a_t m = 2E_1\approx0.52$. We cannot project it out with $\sigma_{00}$
(and its Lanczos time-shifts) because the free triangle $C_{12}=\langle\sigma_{00}\sigma_{00}^2\rangle=0$
by $\sigma_3$-hermiticity ($\pm\lambda$ symmetry). A **local** furnished scalar is no good either:
$\tau(x,x)=\tfrac12 I$ exactly (contact), so $\bar\psi\tau(x,x)\psi=\tfrac12\sigma_{00}$ (verified).

We need an operator that (i) puts $\bar\psi,\psi$ in **different** $(\ell,m)$ modes so it overlaps
$(2,2)$, and (ii) has **nonzero mixing** with $\sigma_{00}^2$. NM's proposal: the covariant
**point-split** scalar on the triangles (one-link separation):
$$
O_\text{ps} = \sum_{\langle i j\rangle} A_{ij}\; \bar\psi_i \,\Omega_{ij}\,U_{ij}\,S\,\psi_j ,
$$
$\Omega_{ij}$ = spin connection (parallel transport of the spinor $j\to i$), $U_{ij}$ = U(1) gauge link
(=1 free), $S$ = scalar kernel (=1 for PS). This is the hopping structure of $D_W$ (`dirac_simp.h`):
$\Omega(i,j)=\cos\tfrac{\omega_{ij}}{2}\sigma_0 - i\sin\tfrac{\omega_{ij}}{2}\sigma_3$.

## Files

- `geom_hopping_claude.py` (NEW): load `omega_n{L}.dat`, `alpha_n{L}.dat`, build the neighbor list,
  $\Omega_{ij}$ (and $\gamma_{ij}=\cos\alpha\,\sigma_1+\sin\alpha\,\sigma_2$), and the point-split
  vertex $W^\text{ps}$ (a $2N_s\times2N_s$ spatial spin matrix).
- `point_split_test_claude.py` (NEW): build $\Phi^\text{ps}[t]=V(t)^\dagger W^\text{ps}V(t)$; test
  (a) mixing $\langle O_\text{ps}\,\sigma_{00}^2\rangle$ (must be $\neq0$), (b) 2-point
  $\langle O_\text{ps}O_\text{ps}\rangle$ effmass (what states, is $(2,2)$ there).
- `two_meson_gevp_pointsplit_free_claude.py` (NEW): $\{1,\sigma_{00},O_\text{ps},\sigma_{00}^2\}$ GEVP.

## Chunks

1. **Geometry loader + vertex** (`geom_hopping_claude.py`). Load omega/alpha, `nns`, build
   $\Omega_{ij}$, assemble $W^\text{ps}[2N_s,2N_s]$ with 2x2 blocks $A_{ij}\Omega_{ij}$ on links.
   Weight $A_{ij}$: start uniform ($Y_{00}$ per directed edge, both directions -> Hermitian-ish);
   refine to a triangle-area measure if needed. Self-test: $\Omega$ unitary, $W^\text{ps}$ shape.
2. **Operator test** (`point_split_test_claude.py`). $\Phi^\text{ps}[t]$; mixing triangle vs $\sigma^2$;
   2-point effmass. Decision gate: mixing nonzero AND 2-pt shows a state near $2E_1\approx0.52$.
3. **GEVP** (`two_meson_gevp_pointsplit_free_claude.py`): $\{1,\sigma_{00},O_\text{ps},\sigma_{00}^2\}$,
   identity carries vacuum, full $G_4$ (both terms) for $\sigma^2$. Target: a level near $2m_\sigma=0.756$
   (two-meson) once $(2,2)$ is projected onto $O_\text{ps}$.

## Open questions
- Link measure $A_{ij}$ for a clean $\ell=0$ projection (uniform vs dual-triangle-area). Test both.
- Whether $S=1$ (pure $\Omega$ transport) suffices, or the Wilson $(-r\sigma_0+\gamma)$ structure is
  wanted (that would make it literally $\bar\psi D_W^\text{hop}\psi$). Start $S=1$.
- Hermiticity: sum over both edge directions, or symmetrize the correlator matrix (solve_gevp does the
  latter already).

## OUTCOME (2026-09-08)

- **Point-split $O_\text{ps}=\bar\psi\Omega U\psi$ FAILED** (both $\Omega$-only and Wilson $-r\sigma_0+\gamma$):
  mixing $\langle O_\text{ps}\sigma^2\rangle\approx0$ (noise), 2-pt $\to0.38$ (ground $(1,1)$). Reason: any
  $\sigma_3$-even geometric vertex has $\mathrm{Tr}[\Gamma GGG]=0$ ($\sigma_3 G\sigma_3=-G$), so the
  $\sigma$-$\sigma^2$ triangle still cancels. `geom_hopping_claude.py` is validated
  ($\Omega$ unitary, $W^\text{ps}$ Hermitian) but the operator doesn't help.
- **$O_A=\bar\psi\tilde\tau\psi$ SUCCEEDED** (mode vertex $\Phi_A=\Phi_{00}\tilde\tau$): mixing
  $\langle O_A\sigma^2\rangle\neq0$ ($\to0.56$), 2-pt $\to0.56$ ($(2,2)$), $\langle O_A\sigma_{00}\rangle=0$
  (orthogonal). The propagator breaks $\sigma_3$-hermiticity ($\sigma_3\tilde\tau\sigma_3=-\tilde\tau-1$).
- **GEVP $\{1,\sigma_{00},\sigma_{00}^2,O_A\}$ (`two_meson_gevp_OA_free_claude.py`) RESOLVES the two-meson:**
  level0=vac, level1=$m_\sigma$=0.378, level2=$(2,2)\approx0.56$, **level3=$2m_\sigma$=0.756 clean plateau.**
- NEXT: L2 confirmation; then the interacting ensembles + couple $F^2$ (the CP $0^{++}$ mixing).

## Refs
- `qed3int_v3-4.pdf` Eq (5.1)-(5.6) (two-term $\sigma_{PS/FS}$); `qed3_v2-6.pdf` App C (free spectrum,
  Eq C.17 $\lambda=n+|m|+1/2$); `dirac_simp.h` (hopping, $\Omega$, $\gamma$, kappa).
