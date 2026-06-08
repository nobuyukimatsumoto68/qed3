# Axial ylm tower -- implementation plan

## Goal

Add the **axial** spherical-harmonic tower $C_{A+-,\ell}(t_0\to t)$ ($\ell=0,1,2$, m-summed), mirroring the
existing VECTOR ylm tower, to the unified correlator program.  Output key
`h{h}/t0_{b}/axial/ylm/Apm/l{l}` (length $N_t$, $\Delta t$-indexed), alongside the existing
`axial/{tp,sp}/Apm`.  (Currently ylm is vector-only; axial is tp/sp only.)

## Physics (mirror vector ylm; ride the axial tp sink)

The axial $C_{A+-}$ legs (PDF Eqs. 3.55/3.68-3.69; valence cases of Sec. 1.1): forward leg
$\phi'_A=X^{-1}\eta$ ($X=D_\text{ov}$ flavor, $D_m$ else), GW factor $\chi=(1-D_\text{ov})\phi'_A$, sink
applies $K^\dagger$.  The m-summed weight is the SAME $W_\ell(n)=A_n/\kappa_t\sum_m Y_{\ell m}$ as the vector
tower.  So the axial ylm tower is:

$$
\psi^A_{\ell}(t_0) = X_\text{sink}^{-1}\,(1-D_\text{ov})\sum_n W_\ell(n)\,K^\dagger(n,t_0)\,\eta,\qquad
\Phi^A_\ell(t) = \sum_n W_\ell(n)\,K^\dagger(n,t)\,\chi,
$$
$$
C_{A+-,\ell}(t_0\to t) = \big(\psi^A_\ell(t_0)\big)^\dagger\,\Phi^A_\ell(t).
$$

Note $\sum_n W_\ell(n)(1-D_\text{ov})K^\dagger\eta=(1-D_\text{ov})\sum_n W_\ell(n)K^\dagger\eta=(1-D_\text{ov})\,
\texttt{srcL}[\ell]$, so the source reuses the SAME $\texttt{srcL}[\ell]=\sum_n W_\ell(n)\,\texttt{rho\_tp}[n]$
($\texttt{rho\_tp}=K^\dagger\eta$, cached from the vector $(++)$ source).  The sink tower
$\Phi^A_\ell$ RIDES the axial tp sink loop (which already computes $K^\dagger(n,t)\chi$ for all $t$) at ZERO
extra $K$ applies -- exactly as the vector ylm tower rides the vector tp sink.

Cost added per hit: $N_\ell\times n_{t0}$ extra source solves (small), zero extra $K$ applies for the sink.

## NOT folded

$C_{A+-}$ is not self-reflection-even (it reflects to $C_{A-+}$, Eq. 3.57), so the axial ylm tower is a
single complex channel, plotted unfolded (like axial tp/sp).

## Files

- EDIT `jj_corr_block_t_claude.cu` (the deployed CORR binary): extend the **axial tp** block (source +
  sink) to build the axial ylm tower and write `axial/ylm/Apm/l{l}`.  Reuse buffers `psi_yl`, `PhiLt`,
  `srcL` (free after the vector ylm block).  New accumulator `Ayl` (N_ELL*n_t0 x Nt).
- (Reference `jj_corr_claude.cu` would need the same edit if used; deferred unless requested.)
- EDIT `jj_corr_massless_claude.ipynb` (+ the mF/mP clones, and the notebook's axial block): loader
  `load_conn('axial/ylm/Apm/l{l}')` -> a per-$\ell$ axial tower plot (unfolded, signed).

## Ordered chunks

1. **Program.** Axial tp source: add the $\ell$ loop building `psi_yl[IYL(l,b)]` $= X_\text{sink}^{-1}
   (1-D_\text{ov})\,\texttt{srcL}[\ell]$ (single solves; flavor/parity/massless op set as axial tp).
   Axial tp sink: zero `PhiLt`, accumulate `PhiLt[IPL(l,t)] += W_ell[l][n]*kphi` in the n-loop, then after
   the n-loop `Ayl[IYL(l,b)][dt] += psi_yl[IYL(l,b)].dag(PhiLt[IPL(l,t)])`; write `axial/ylm/Apm/l{l}`.
   Files: `jj_corr_block_t_claude.cu`.  (User compiles/validates.)
2. **Notebook.** Axial ylm loader + per-$\ell$ unfolded signed plot in the massless notebook (and clones).

## Validation

Free massless: the axial ylm tower should fall like the vector one in the conformal regime (per-$\ell$
decay), single complex channel; $\ell=0$ small.  Cross-check `axial/tp/Apm` unchanged (bit-identical to the
pre-edit run, since the axial tp source/sink lines are untouched -- only ADD the ylm accumulation).
