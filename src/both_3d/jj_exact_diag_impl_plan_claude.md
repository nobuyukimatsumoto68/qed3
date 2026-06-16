# Impl plan — insertion-diagonal deterministic EXACT current (jj_exact_diag_deter_claude.cu)

Goal: a deterministic exact conserved-current (K) jj that is **insertion-diagonal** (matches disp/loc and
the production stochastic `jj_corr_block_t`, and reproduces Eqs. 4.31/4.28), unlike the existing
`jj_kbuild_exact`+`jj_contract_deter` which contract the SUMMED $K^P$ -> a q=0 DOUBLE sum.

Estimator (vector tp + sp), per qed3int_v2-12 (4.31 tp same-site / 4.28 sp via triangle id):
$$
C^P(t_0\!\to\! t)=\sum_a w^P_a\,\mathrm{tr}\big[K(a,t_0)\,P\,K(a,t)\,P\big],\qquad
\text{tp: }a\in\text{sites},\ w_{tp}=\tfrac{A_n}{\kappa_t^2};\quad
\text{sp: }a\in\text{links},\ w_{sp}=\tfrac{\ell\ell^*}{\kappa^2}.
$$
$K(a,t)$ = full overlap conserved current (ConservedCurrent `kop` + `op_K`, multishift resolvent) =
the SAME kernel/weights as `jj_kbuild_exact` -- only the contraction is diagonal.

Implementation (merge kbuild GPU work + contract time-shift):
- Load $P=D_m^{-1}$ (`Dm_inv`).
- For each proj, each insertion $a$: `kop.set_temporal/spatial(U,0,a)`; build $B_a=K(a,0)\,P$ by
  applying `op_K` to each COLUMN of $P$ (N applies; == kbuild cost $N\times n_\text{ins}$). One $B_a$
  ($N\times N$) live at a time.
- Free field => $C^P(dt)=\sum_a w_a\,\texttt{conn\_shift}(B_a,dt)$ (antiperiodic time-shift, reused from
  jj_contract_deter); disc $=\sum_a w_a\,\mathrm{tr}(B_a)$.
- Output `data_<ESNID>/corr_deter_exactdiag_L<L>/corr.<k>.h5`, keys `h0/t0_b/{tp,sp}/{Vpp,Vmm}` +
  `h0/disc/{tp,sp}/J` (jj-style; DISTINCT dir from the double-sum `corr_deter_exact_L*`).

CLI mirrors jj_kbuild_exact (mass via --mass-re/im for the P dir; --prop-file/--out-tag; --n-t0).
L=1,2 feasible (same op_K count as kbuild; memory one $N\times N$).

Run: `tmp_exactdiag_claude.sh` (lattice prop, L=1,2).

Notebook `comp_trio_jj_claude.ipynb`: tp = loc `s3` vs exactdiag `tp`; sp = loc `s1`+`s2` vs
disp `sp` vs exactdiag `sp`. Renormalize to a common point (match at one separation); the headline is
the **sign** and shape of exact sp vs disp/loc.

## Addendum (2026-06-11) — `--sum` mode: diagonally-n-summed exact (Eq.4.29)

Fold a `--sum` flag into `jj_exact_diag_deter_free_claude.cu` (NOT a new file).  When set, compute the
area-weighted INSERTION-DIAGONAL sum (Eq.4.29), free + deterministic, L=1,2:
$$
G^{tp}=\frac{1}{4\pi}\sum_{n\in\text{sites}} A_n\,\mathrm{tr}[K_{ov\kappa,n}(0)PK_{ov\kappa,n}(t)P],\quad
G^{sp}=\frac{1}{4\pi}\sum_{\ell\in\text{links}} A_\ell\,\mathrm{tr}[K_{ov\kappa,\ell}(0)PK_{ov\kappa,\ell}(t)P]
$$
with $A_n=$`dual_areas`, $A_\ell=$`link_volume`, $K_{ov\kappa}=K/\kappa$ (`insertion_kappa`), $1/(4\pi)$ folded by
`write_corr`.  Per insertion: build dense $K$ (N `op_K` applies) /\kappa, $A=K\cdot P$ (`matmul_A`),
accumulate $A_n\,$`conn_shift`$(A,dt)$.  BUILD-USE-DISCARD (no cache: all-n caching is ~300 GB at L=2).
Output `corr_deter_exactsum[_<tag>]_L<L>`, keys `h0/t0_b/{tp,sp}/{Vpp,Vmm}` + disc.  Cost = (n_sites+
n_links)*N solves: L=1 ~min, L=2 ~hours (flagged).  CHECK: ratio $G^{sp}/G^{tp}\to-2$ (geometric sum
restores the D-1 directions; single-insertion exact1 was -1).  Single-insertion path (default, no --sum)
is unchanged.
