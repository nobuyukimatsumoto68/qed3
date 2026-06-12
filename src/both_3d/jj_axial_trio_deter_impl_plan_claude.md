# Axial trio (deterministic, free-field) -- implementation plan

Goal: reproduce `comp_trio_jj_claude.ipynb` (loc vs disp vs exact, diagonally-$n$-summed, normalized
currents, lattice/continuum prop, $L=1,2$) but for the **AXIAL** current $C_{A+-}$, instead of the
vector $C_{V++}$.

Reference: `qed3int_v2-13.pdf` (src `qed3int_v2(6)/main.tex`), Sec. 3.3.2 massless.  Algorithm provenance
of the dense free-field contraction: our own `jj_exact_diag_deter_free_claude.cu` (`conn_shift`,
antiperiodic time-shift); axial operator structure validated (stochastically) in
`jj_corr_block_t_claude.cu` (audited in `jj_corr_pdf_audit_claude.md`).

## Physics

### Exact axial -- Eq.(3.55)
$$
C_{A+-,c}(t_0\to t)=\mathrm{tr}\!\big[K(t_0)\,(1-D_\text{ov}^\dagger)\,D_\text{ov}^{\dagger-1}\,
K^\dagger(t)\,(1-D_\text{ov})\,D_\text{ov}^{-1}\big].
$$
Vector $C_{V++}$ (Eq.3.52) is the same with no GW factors and $K(t)$ (not $K^\dagger$) on the sink and
$D^{-1}$ (not $D^{\dagger-1}$) on the $t_0$ leg:
$$
C_{V++,c}(t_0\to t)=\mathrm{tr}\!\big[K(t_0)\,D_\text{ov}^{-1}\,K(t)\,D_\text{ov}^{-1}\big].
$$
So the axial differs from the vector in three places: (i) GW factors $(1-D_\text{ov}^\dagger)$ on $t_0$,
$(1-D_\text{ov})$ on $t$; (ii) the $t_0$-leg propagator is $D^{\dagger-1}=P^\dagger$ (not $P$); (iii) the
sink kernel is $K^\dagger(t)$ (not $K(t)$).

### Local / displaced axial -- "replace $K\to\sigma$ (or $W_d$), rip off $(1-D)$"
Drop both GW factors (they are $O(a)$ in the local/displaced reductions) and replace the conserved
kernel by the bare insertion:
$$
C^{\rm loc}_{A+-}=\mathrm{tr}\!\big[\sigma^a(t_0)\,P^\dagger\,\sigma^a(t)\,P\big],\qquad
C^{\rm disp}_{A+-}=\mathrm{tr}\!\big[W_d(t_0)\,P^\dagger\,W_d^\dagger(t)\,P\big].
$$
($\sigma^a$ is Hermitian so $\sigma^{a\dagger}=\sigma^a$; $W_d$ is the link kernel $-P_{xy}U_{xy}$ and is NOT
Hermitian, so the sink carries the genuine dagger $W_d^\dagger$.)  Compared to the **vector** loc/disp
(`jj_local_deter`, `jj_disp_deter`), the ONLY changes are: the $t_0$-leg propagator $P\to P^\dagger$ (i.e.
$D^{\dagger-1}$), and the sink kernel daggered.  $P^\dagger_{ab}=\overline{P_{ba}}$ -- no extra solve, a
conjugate-transpose lookup of the loaded dense $P=D_m^{-1}$.

### Reflection / reality
$C_{A+-}$ is NOT self-reflection-even: it reflects to $C_{A-+}(t_0\to t)=C_{A+-}(t\to t_0)$ (Eq.3.57).
Hence it is a SINGLE complex channel (no $V$mm$=$conj partner).  Plot $\mathrm{Re}$ as the signal, carry
$\mathrm{Im}$ as a diagnostic (unfolded, like the stochastic `axial/{tp,sp}/Apm`).

## Dense free-field implementation (exact)

Established framework (`jj_exact_diag_deter_free`): build dense $K(a,0)$ col-by-col ($N$ `op_K` applies),
$A_0=K\,P$, and `conn_shift(A0,dt)`$=\mathrm{tr}[A_0\,S_{dt}A_0 S_{-dt}]=\mathrm{tr}[K(0)PK(dt)P]$ using the
ANTIPERIODIC time shift $K(dt)=S_{dt}K(0)S_{-dt}$ and $S_{dt}P S_{-dt}=P$.

Axial generalization (two distinct matrices):
- $G\equiv(1-D_\text{ov})$ dense, built once via `op_oneMinusD` col-by-col ($N$ applies, mass- and
  insertion-INDEPENDENT -> cache `data_free_Gcache_L<L>/G.h5`, reuse across insertions and tp/sp).
- $K^\dagger(a,0)$ dense $=$ conjugate-transpose of the cached dense $K(a,0)$ (no extra applies).
- Source-side: $B_\text{src}=K\,G^\dagger\,P^\dagger$ (two matmuls; $G^\dagger,P^\dagger$ are
  conj-transposes).
- Sink-side: $B_\text{snk}=K^\dagger\,G\,P$ (two matmuls).
- $C_{A+-}(dt)=$ `conn_shift2(B_src,B_snk,dt)`$=\mathrm{tr}[B_\text{src}\,S_{dt}B_\text{snk}S_{-dt}]
  =\mathrm{tr}[K(0)G^\dagger P^\dagger\,K^\dagger(dt)G P]$, which is Eq.(3.55) with $t_0=0$.

`conn_shift2(X,Y,dt)` is the existing `conn_shift` with the second matrix $Y$ instead of $X$:
`s += X[p,q]*(s_p s_q)*Y[shift(q,-dt),shift(p,-dt)]`.  Validity of the whole-matrix shift on $B_\text{snk}$:
$B_\text{snk}(t)=K^\dagger(t)GP=S_t(K^\dagger(0)GP)S_{-t}$ since $G,P$ are time-translation invariant.

Cost (per insertion): reuse the K cache (0 op_K applies on hit) + build $G$ once + 4 dense $N\times N$
matmuls.  $L=1$ ($N{=}3072$) trivial; $L=2$ ($N{=}10752$) minutes for the full `--sum`.  $L=4$ ($N{=}41472$,
$N^2{\sim}1.7$e9) NOT targeted here (trio is $L=1,2$).

## Local / displaced implementation
`tr_WPWP` (O(16) loc / O(64) disp dense-$P$ lookups) generalized so the FIRST propagator uses $P^\dagger$:
```
tr_WPWP_axial(E0, Et, P):   // Et = W^dag at sink (loc: sigma unchanged; disp: daggered entries)
  s += e0.v * conj(P[et.i*N + e0.j]) * et.v * P[et.j*N + e0.i]   // Pdag[e0.j,et.i] = conj(P[et.i,e0.j])
```
- loc: `Et = local_W_sigma(...)` unchanged (sigma Hermitian); only the first-propagator lookup conjugated.
- disp: build $W_d^\dagger$ at the sink (swap $i\leftrightarrow j$, conjugate $v$ of `build_W_disp`).
Channels, weights, ylm tower, `--sum`/`--ins`, output layout: IDENTICAL to the vector codes; only the
dataset name `Vpp`->`Apm` (single complex channel; no `Vmm`).

## Files (proposed: NEW side-by-side copies; vector codes untouched)
- `jj_exact_axial_deter_free_claude.cu`  (copy of `jj_exact_diag_deter_free_claude.cu`)
  -> `corr_deter_exactsum_axial[_<tag>]_L<L>`, keys `h0/t0_b/{tp,sp}/Apm`, `.../ylm/l{l}/Apm`.
- `jj_local_axial_deter_claude.cu`       (copy of `jj_local_deter_claude.cu`)
  -> `corr_deter_local_axial[_<tag>]_L<L>`, keys `h0/t0_b/{s1,s2,s3}/Apm`, `.../ylm/l{l}/Apm`.
- `jj_disp_axial_deter_claude.cu`        (copy of `jj_disp_deter_claude.cu`)
  -> `corr_deter_disp_axial[_<tag>]_L<L>`, keys `h0/t0_b/{tp,sp}/Apm`, `.../ylm/l{l}/Apm`.
- `comp_trio_jj_axial_claude.ipynb`      (copy of `comp_trio_jj_claude.ipynb`; loaders -> `_axial` dirs,
  `Vpp`->`Apm`; same tp/sp/ratio/ylm cells; add an Im-diagnostic note).
- run scripts: copies of the existing `tmp_exactsum_free_claude.sh` / `tmp_summed_fast_claude.sh` etc.

## Ordered chunks
1. **Exact axial** (`jj_exact_axial_deter_free_claude.cu`): `conn_shift2`, $G$ cache, $B_\text{src}/B_\text{snk}$,
   `--sum` tp/sp + ylm, write `Apm`.  (single-`--ins` path too, for spot checks.)  User compiles/runs $L=1$.
2. **Local axial** (`jj_local_axial_deter_claude.cu`): `tr_WPWP_axial`, $P^\dagger$ first leg, write `Apm`.
3. **Disp axial** (`jj_disp_axial_deter_claude.cu`): `tr_WPWP_axial` + $W_d^\dagger$ sink, write `Apm`.
4. **Notebook** `comp_trio_jj_axial_claude.ipynb`.

## Resolved decisions (2026-06-12) + status: IMPLEMENTED (awaiting user compile/run)
1. File strategy: NEW side-by-side `_axial` copies (vector codes byte-identical).  DONE.
2. Scope/order: all three programs in one pass, then the notebook.  DONE.
3. Channel comparison: mirror the vector trio (plot $\mathrm{Re}\,C_{A+-}$, ratio $G^s_A/G^t_A$, ylm rates);
   $\mathrm{Im}$ carried as a diagnostic cell.  $C_{A-+}$ reflection NOT computed (reconstruct downstream
   if needed).  DONE.
4. CFT anchor: free-theory axial $=$ vector -> the notebook draws the same $-1$ / $e^{-2t},e^{-3t}$ / ratio
   $2.4$ references; the data reveals whether the axial follows them.

Files written: `jj_exact_axial_deter_free_claude.cu`, `jj_local_axial_deter_claude.cu`,
`jj_disp_axial_deter_claude.cu`, `comp_trio_jj_axial_claude.ipynb`, run scripts `tmp_axial_fast_claude.sh`
(loc+disp, GPU0, L=1,2) and `tmp_axial_exactsum_claude.sh` (exact `--sum`, GPU1, L=1; add L=2 when ready).
Bug fixed pre-handoff: hoisted `D.update(U)` in the exact program so the $G=(1-D_\text{ov})$ build is valid
on a K-cache-HIT path.  Memory note: exact `--sum` at $L=2$ holds ~33 GB host (9 ylm Sigma + G/Gdag/P/Pdag +
per-insertion working) -- same order as the vector exactsum; $L=1$ is trivial (run that first).

## Validation expectations
- `Re` tp/sp/ratio: if axial follows free vector, ratio $G^s_A/G^t_A\to-1$, tp$<0$, sp$>0$; loc/disp/exact
  agree after the $dt$=5 match.  ylm rates $\to(2,3)$, $g_2e^{3t}/g_1e^{2t}\to2.4$, $g_0\to0$.
- `Im` diagnostic: free massless even-$N_f$ -> consistent with noise/cutoff (jj_reality_tsymmetry_claude.md).
  A stable nonzero `Im` would flag a genuine axial phase or a source/sink ($K$ vs $K^\dagger$) asymmetry.
