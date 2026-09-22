# Findings: the $\Omega$-rotated free-kernel construction of the $(2,2)$ interpolator is obstructed for $U(1)/S^2$

Author: session "Explore: $\Omega$-rotated free approx" (qed3-79). Date: 2026-09-18.
Collaborators: "Fin: Two-meson" (handoff, perambulator/GEVP plumbing, ckpoint reader), "main-thread [5c528b]" (SU(N) flow + Landau-fix method author).
Companion handoff: `omega_rotated_free_kernel_handoff_claude.md`.

## Verdict (one line)

Dressing a free $\ell=3/2$ shell projector with a Landau gauge rotation, $K=\Omega^\dagger P^{\text{free}}\Omega$, **cannot restore the interacting $(2,2)$ state**. The mechanism is fundamentally obstructed -- confirmed independently by the $U(1)$ Hodge argument here and by [5c528b]'s "topological floor" in $SU(N)$. Use Route 1 (`sigma2_shell22_gevp_claude.py`, interacting L2 $(2,2)\approx0.65$) for the physics.

## What was tried

Target: a gauge-invariant single-meson interpolator overlapping the $0^{++}$, $\ell=0$, $\Delta=4$ $(2,2)$ state (both constituents in the $\ell=3/2$ spatial shell on the $S^2$ timeslice). Construction (see handoff):
$$
K_{(2,2)}(x,y) = \Omega^\dagger(x)\,P_{3/2}^{\text{free}}(x,y)\,\Omega(y),
$$
with $\Omega=e^{i\phi}$ the per-timeslice Landau fixing rotation of a (moderately flowed) config, $P_{3/2}^{\text{free}}$ the exact free shell projector. Correlator via the distillation machinery: $\Phi_K(t)=V(t)^\dagger K(t)V(t)$, then $C(t,s)=-\mathrm{Tr}[\Phi_K(t)\tau(t,s)\Phi_K(s)\tau(s,t)]$.

Convention (pinned with Fin + [5c528b]): $\Omega$ = fixing rotation, smooth $=\Omega U\Omega^\dagger$, $\psi\to\Omega\psi$; the sandwich `phi=Omega*in; Pfree(phi); adj(Omega)*out` realizes $\Omega^\dagger P^{\text{free}}\Omega$. Sign pinned by the convention-free rule: after building $\Omega$, the Landau functional must drop.

## The obstruction (rigorous)

For $U(1)$, $\Omega=e^{i\phi}$ with $\phi$ solving the Poisson equation
$$
\text{Lap}\,\phi = -\,\text{div}\,a,
$$
so $\phi$ (hence $\Omega$) depends ONLY on $\text{div}\,a$, the longitudinal (pure-gauge) part of the connection $a_i=\arg U_i$.

The spatial Wilson flow force is the co-derivative of the curvature, $\dot a = \delta F$ with $F=da$. Then
$$
\text{div}(\dot a) = \delta\,\delta F = \delta^2 F = 0,
$$
so the flow preserves $\text{div}\,a$ EXACTLY. Verified numerically: max$|\text{div}(a_\text{raw})-\text{div}(a_\text{flowed})|=$ 1e-15 per slice while the plaquette action drops $S:31\to0.09$ (flow $\tau=3$).

Therefore $\Omega$ is flow-invariant: for $U(1)$, "flow then Landau" $\equiv$ "Landau alone." A site-diagonal phase $e^{i\phi}$ removes only the longitudinal/pure-gauge distortion of the free projector; the TRANSVERSE (physical) fluctuations -- which are what actually distort the interacting shell -- are untouched by any gauge rotation. At $g^2=1$ this transverse roughness is $O(1)$ physical curvature (not a UV artifact), so the free shell simply does not map to the interacting shell under $\Omega$.

[5c528b] confirms this is the same "topological floor" they hit in $SU(N)$: a gauge rotation (abelian or not) provably cannot carry gauge-invariant/physical content. The flow's real value in their construction was transverse smoothing, not improving $\Omega$; the "flow moves the Landau frame" effect (present in $SU(N)$ via non-commutativity, absent in $U(1)$ by Hodge) was incidental.

## Empirical results (interacting L2 $g^2=1$, $N_f=2$, truncated $N_v=24$ basis)

Reference numbers: $m_{PS}=0.353$, $2m_{PS}=0.705$, Route-1 shell-GEVP $(2,2)\approx0.65$, free-L2 $2E_{3/2}\approx0.676$.

- Naive $\Omega=1$ insertion $V^\dagger P^{\text{free}}V$: effmass falls monotonically to $\approx0.37$ ($=m_{PS}$) by $dt\sim12$ -- the free projector loses shell character on interacting configs and picks up the ground state.
- $\Omega$-dressed $K=\Omega^\dagger P^{\text{free}}\Omega$: effmass $\approx0.37$, BYTE-IDENTICAL to naive, and identical at flow $\tau=0,0.2,1,3$ (because $\Omega$ is flow-invariant).

So the dressing does essentially nothing; both leak to $m_{PS}$.

## Validated, reusable assets (all in `src/production`)

- Position-space free projector $P^{\text{free}}$: built as $V_\text{free}\,\Phi_\text{shell}\,V_\text{free}^\dagger$ from the complete free basis (`data_free/distill_Nv84`, $N_v=2N_s=84$). Exact rank-8 projector (trace 8.000, hermitian & idempotent to 1e-16, eigenvalues $\{1^{\times8},0^{\times76}\}$), t-independent to 6e-16. Mode-space version reproduces free L1 $2E_{3/2}=0.556$. The $\ell=3/2$ window is ranks $[4,12)$ at BOTH L1 and L2 (free $|\text{eig}(M)|$ clusters $4,8,12$ do not smear).
- Tested ckpoint reader `ckpoint_reader_claude.py` (from Fin): spatial link phases, edge map (`links_n<L>.dat`, NOT `omega_n<L>.dat`), graph Laplacian, divergence, face flux. Edge order validated by mean plaquette $\langle\cos\rangle=0.664>0$.
- Spatial Wilson flow on the $S^2$ triangulation (gradient descent on $\sum_\text{faces}(1-\cos\text{flux})$; each link in exactly 2 faces): $S:31\to0.09$ at $\tau=3$, flux sectors preserved.
- $U(1)$ Landau/Poisson $\theta$-solve $\text{Lap}\,\phi=-\text{div}\,a$; sign +1 drops the Landau functional on every slice.
- Branch note: at $g^2=1$ the raw spatial phases have $|a_i|<\pi$ on every link, so the arg-branch is unambiguous before any flow (re-check at stronger coupling).

## Salvage options and recommendation

1. Route 1 (`sigma2_shell22_gevp_claude.py`): shell operators from the $M=\tau(t,t)-1/2$ eigenspaces $[0,4],[4,12],[12,24]$ + a GEVP that orthogonalizes $m_{PS}$ out $\to$ interacting L2 $(2,2)\approx0.65$. Physics goal ALREADY met.
2. Flowed-link projector $P[D_W[U_\text{flow}]]$: put the flow's transverse smoothing into the FERMION object. Gauge-invariant, sounder shell operator. Pitfalls (Fin + [5c528b]): converges back onto Route 1 (still needs the GEVP as a diagonal); re-enters a lattice eigenbasis unless the flowed shell overlaps the free-continuum shell well enough to apply $P^{\text{free}}$ without diagonalizing (test that overlap explicitly); flow-time window can wash out the $(2,2)$; $g^2=1$ roughness is $O(1)$ physical so expect only partial help. Cross-check value, not new capability.
3. [5c528b]'s alternative: build $(2,2)$ directly from covariantly-smeared / covariant-derivative fermion fields on the (flowed) links -- angular structure AND gauge invariance without asking a free projector to carry interacting physics. The standard smeared/derivative-source route; genuinely different, more work.

Recommendation: Route 1 for the physics (done). If an independent handle is wanted, prefer option 3 over option 2 (option 2 mostly reproduces Route 1).

## References

- Chester, Pufu, arXiv:1603.05582 (the $0^{++}$ $(\bar\psi\psi)^2$ vs $F^2$ physics context).
- Fourier-accelerated Landau gauge fixing: Davies et al., PRD 37 (1988) 1581.
