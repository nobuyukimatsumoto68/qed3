# Implementation plan — check Eq. (4.31) for $G_t$ via the continuum propagator

Source of equations: `jj_cft_check_equations_claude.md` (records qed3int_v2-11 Eqs. 1.24, 3.6-3.11,
4.26-4.31). Algorithm/method: standard free-fermion-loop current correlator; current normalization
and the radial-quantization CFT form are from qed3int_v2-11 Sec. 4 (no external algorithm to cite
beyond the project's own derivation).

## Goal

Verify the temporally projected conserved-current two-point function, Eq. (4.31):
$$G_t(t;\hat n,\hat n) = -C_J\,e^{-\Delta t}(1-e^{-t})^{-2\Delta},\qquad \Delta=2,$$
by contracting the (validated) continuum eigenbasis propagator `G_prop`
(`cont_prop_eigbasis_claude.cpp`, Boost $\xi$) into the vector current $j_V^a=\eta^\dagger\sigma^a\xi-\xi^\dagger\sigma^a\eta$:
$$f_{\rm cyl}^{ab}(t;\hat n,\hat n) = -\,\mathrm{tr}[\sigma^a G(x_1,x_2)\sigma^b G(x_2,x_1)],\quad
  G_t=f^{33},\ \ G_s=f^{11}+f^{22}.$$
Checks: (i) $t$-shape $G_t\propto e^{-2t}(1-e^{-t})^{-4}$ (ratio $G_t/[\,e^{-2t}(1-e^{-t})^{-4}]$
constant), (ii) $\hat n$-independence (several $(\theta,\phi)$ give the same $G_t$), (iii)
$G_s/G_t=-2$, (iv) agreement with the closed form $G_t=2c_3(t)^2=\tfrac{1}{8\pi^2}e^{-2t}(1-e^{-t})^{-4}$.

## Files

- NEW `jj_cft_Gt_check_claude.cc` — standalone; copies `xi_mn` (Boost), `c2_mn`, `cnm2_mn`,
  `lambda_mn`, `G_prop` from `cont_prop_eigbasis_claude.cpp`; adds $\sigma^a$, a 2x2 trace
  contraction, the $t$-sweep, and the comparisons. Pure host (GSL + Boost). NO CUDA.
- NEW `run_jj_cft_Gt_check_claude.sh` — `g++ -x c++` build (+ Boost `-I`), run, tee to log.
- Equations: `jj_cft_check_equations_claude.md` (done).

## Chunks

1. **Propagator plumbing.** Copy `xi_mn`(Boost)/`c2_mn`/`cnm2_mn`/`lambda_mn`/`G_prop` verbatim.
   Add `sigma[3][2][2]` (Pauli) and static helpers `mat2_mul`, `mat2_trace`.
   Files: `jj_cft_Gt_check_claude.cc`.

2. **Contraction.** `fcyl_ab(theta,phi,t,nmax) -> cd F[3][3]`:
   `G12 = G_prop(theta,phi,t; theta,phi,0)`, `G21 = G_prop(theta,phi,0; theta,phi,t)`;
   `F[a][b] = -trace(sigma[a]*G12*sigma[b]*G21)`. (Diagonal expected: $F=\mathrm{diag}(-2,-2,+2)c_3^2$.)
   Files: `jj_cft_Gt_check_claude.cc`.

3. **Sweep + compare.** For $t$ in a moderate grid (e.g. 0.4..4.0), at a reference $\hat n$:
   print $G_t=F^{33}$, $G_s=F^{11}+F^{22}$, $G_s/G_t$, $\mathrm{ref}=e^{-2t}(1-e^{-t})^{-4}$,
   $G_t/\mathrm{ref}$ (should be constant $=1/(8\pi^2)$), and the analytic $2c_3(t)^2$.
   Also evaluate $G_t$ at 3 different $(\theta,\phi)$ to confirm $\hat n$-independence.
   Verdict: ratio-constancy (rel. spread) PASS, $G_s/G_t\to-2$ PASS, $\hat n$-spread PASS.
   Files: `jj_cft_Gt_check_claude.cc`.

4. **(Deferred)** general $n_1\ne n_2$ via Eq. (4.26) and the off-site `G_prop` blocks (the angular
   $\Lambda$ structure). Not needed for (4.31); sets up the later $G_s$/$Y_{\ell m}$-tower checks.

## Numerics

- `nmax`: Boost $\xi$ is stable, so use `nmax=100`. Convergence of $c_3(t)=\frac1{4\pi}\sum(n+1)e^{-(n+1)t}$
  needs $(n_\text{max}+1)t\gg1$; at $t=0.4$, $e^{-101\cdot0.4}=e^{-40}$ — fine. Avoid $t\to0$ (the
  $(1-e^{-t})^{-4}$ UV blow-up); grid starts at $t\approx0.4$.
- Single thread; `g++ -x c++ -I../../qfe_mod/include/boost_1_86_0 ... -lgsl -lgslcblas -lm`.

## Open questions / assumptions to confirm (flagged in the .cc)

1. Temporal projection $e^a_3 e^b_3 \to \sigma^3$ (the $a=b=3$ component), justified by the same-site
   block being pure $\sigma_3$ in this frame. If your $e^a_3$ convention differs, the mapping changes.
2. Current = conserved **vector** $j_V$; overall constant $C_J$ (incl. sign) absorbed — we check
   shape / $\hat n$-independence / $G_s/G_t=-2$, not the absolute $C_J$.
3. $t$ (cylinder time) $=$ `G_prop`'s $\tau$; $\Delta=2$.
