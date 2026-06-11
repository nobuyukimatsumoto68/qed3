# Impl plan: update the L=1 free-correlator notebook to the $C_c$ convention

## Goal / convention

Bring `jj_free_correlators_claude.ipynb` (L=1) in line with the L=2/L=4 notebooks and the agreed sign
convention. The stored/plotted quantity is the **connected** current correlator
$$
C_c^{P}(t_0,t)=\sum_a w^P_a\,\mathrm{tr}\!\big[\,J(a,t_0)\,P\,J(a,t)\,P\,\big]\qquad(P=D_\text{ov}^{-1}),
$$
in the **+ trace convention** (no fermion-loop minus, $1/4\pi$-scaled, insertion-diagonal). The physical
current-current correlator is
$$
\langle JJ\rangle = -\,C_c + C_d ,\qquad C_d=\langle J\rangle\langle J\rangle\ (\approx 0\text{ free}),
$$
so with $\langle JJ\rangle<0$ (tp) and $\langle JJ\rangle>0$ (sp) the expectation for the plotted $C_c$ is
$$
\boxed{\,C_c^{tp}>0,\qquad C_c^{sp}<0\,}\quad(\text{free}, C_d\approx0).
$$
At L=1: **local** satisfies this ($C_c^{tp}>0$, $C_c^{sp}<0$); **exact** has $C_c^{sp}>0$ (wrong sign) —
the open physics question (exact spatial sign), tracked separately, not part of this notebook edit.

## File to modify

- `jj_free_correlators_claude.ipynb` (only this notebook).

## Chunks

1. **Relabel $G\to C_c$ + convention note.** Files: the notebook.
   - Intro markdown (`2b703d8e`): add the $\langle JJ\rangle=-C_c+C_d$ relation and the boxed sign
     expectation; change the CFT line $G_{t,s}\to C_c^{tp},C_c^{sp}$, ratio $C_c^{sp}/C_c^{tp}=-(D-1)=-2$.
   - `plot_reim` (`05bdce71`): titles `Re/Im G_{proj}` -> `Re/Im C_c[{proj}]`.
   - ratio cell (`3933a5b0`): `Gs/Gt` -> `C_c[sp]/C_c[tp]` in label/ylabel/title/print.
   - section markdowns (`e436621f`, `7fddce73`, `36a64d98`): `G`/`|G|`/`$G_s/G_t$` -> `C_c` forms.

2. **(optional) error bars on a stochastic exact curve.** Files: the notebook.
   - Only if we switch/add the stochastic exact (`corr_nt02_nhits64`, per-hit groups `h0..h63`): port the
     unified `load_hits`/`mean_sem`/`errbar_signed` from L2/L4 and draw exact with signed error bars.

## Open question (resolve before editing)

- **Which "exact" at L=1?**  (a) keep the **deterministic basis-loop** $K$ `corr_deter_exact_L1`
  (clean, all-dt, no error bars) — only Chunk 1 needed; (b) switch to the **stochastic**
  `corr_nt02_nhits64` (matches L2/L4, error bars, but noise past dt~15) — Chunks 1+2; (c) **both**
  curves overlaid.
