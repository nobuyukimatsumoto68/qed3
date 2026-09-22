# Magnetic point-to-point GEVP (single fixed point) -- impl plan

qed3-4e / "Fin: axial sp". Resume task from NM: improve the magnetic ($\Psi_{\ell=1}$) extraction via a
GEVP whose partner operator is the **point-to-point current at a single fixed site** $n_0$ (NOT a
sum/split). See `findings_sp_vsh_claude.md` Section 5.

## Physics / goal
The magnetic single-op correlator $C_{00}$ (fully VSH-projected $\Psi_{\ell=1}$) settles late (dt 12-15,
$a_t m\approx0.36$) -- excited-tower contamination. A GEVP needs a second operator with a DIFFERENT
overlap onto the tower. The earlier partners failed because they were isotropic averages/splits:

- time-split $\bar\psi_t\sigma\psi_{t+1}$ and NN-split $\sum_{n'\in NN}\bar\psi(n)\sigma\psi(n')$ both
  DECOUPLED ($C_{11}\to0$): summing over neighbours washes out the directional curl.
- Averaging a "point" partner over all 12 sites reduces it back to the projected op (multiplicity=1).

NM's fix: **pick a single point.** The partner is the magnetic current localized at ONE fixed vertex
$n_0$. A point operator is NOT projected onto a single VSH mode -- it overlaps the FULL tower (all
$\ell$), so its excited-state content differs from the pure-$\ell{=}1$ projected op, and the $2\times2$
GEVP can separate the magnetic ground from contamination. Crucially $C_{11}=M_{n_0 n_0}$ is a healthy
$O(1)$ local magnetic density (real magnetic overlap), so it does NOT decouple like the isotropic splits.

## Operators (magnetic tangent weight $W^a_m(n)$, $a\in\{1,2\}$ local frame, $m\in\{-1,0,1\}$)
- $O_0$ = projected magnetic VSH: $\sum_n W^a_m(n)\,j^a(n)$ (sum over all 12 sites).
- $O_1$ = point magnetic current at fixed $n_0$: $W^a_m(n_0)\,j^a(n_0)$ (single site, no sum).

## Master object: magnetic-directed site-pair matrix (per config, per dt)
$$M_{n_1 n_2}(t)=\sum_{m}\sum_{a,b\in\{1,2\}} W^a_m(n_1)\,W^b_m(n_2)\,f^{ab}_A(n_1,n_2;t)$$
with the AXIAL contraction $f^{ab}_A=-\mathrm{tr}[\sigma^a G_{21}^\dagger\sigma^b G_{21}]$ (spin-daggered
forward leg), exactly as in `interacting_vsh_L1_claude.py` -- but keep the $12\times12$ structure instead
of summing to a scalar. Everything else derives from $M$ (no recompute):
- $C_{00}=\sum_{n_1 n_2}M_{n_1 n_2}$   (projected-projected = the existing magnetic corr)
- $C_{01}=\sum_{n_1}M_{n_1 n_0}$,  $C_{10}=\sum_{n_2}M_{n_0 n_2}$   (proj $\leftrightarrow$ point)
- $C_{11}=M_{n_0 n_0}$   (point-to-point at $n_0$)
$M$ is symmetrized $0.5(M+M^T)$ so $C_{01}=C_{10}$.

## GEVP
Solve $C(t)v=\lambda C(t_0)v$ (2x2), ground = largest $\lambda$, $a_t m=\log(\lambda(t)/\lambda(t{+}1))$.
Orient so ground is positive; hermitize. $n_0$ = vertex 0 (all 12 equivalent by icosahedral symmetry;
pick one -- averaging would collapse the partner). Jackknife over config bins (bin10, kmin20).

Deliverable check: does the GEVP ground plateau EARLIER and tighter than the single-op dt 12-15?
Also print the bare point-to-point $C_{11}$ effmass and $C_{00}$ effmass for reference. Save the full
$12\times12$ binned $M$ so the full site-pair GEVP (or any other point choice) is a cheap post-process.

## Files
- NEW `mag_point2point_gevp_claude.py` (src/production): builds $M_{n_1 n_2}(t)$ from perams, bins,
  forms the 2x2 {projected, point-$n_0$} GEVP, saves `final/analysis_axial/mag_point2point_gevp_claude.npz`
  + `..._claude.png`. Reuses `fs_gevp_point_claude.make_config` and the axial $f$ contraction; GPOF/log
  effmass inline. ENS via env (default Nf2 g1.0 L1 distill_Nv24).

## Open question (flag, do not block)
Basis assumed = 2-op {full-projected $O_0$, single-point $O_1(n_0)$}. If NM instead wants the full
$12\times12$ site-pair GEVP, it is a one-line change on the saved $M$ (post-process). Confirm after seeing
the 2x2 result.
