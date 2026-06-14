# Impl plan -- spacetime-averaged condensates (Eq. 1.23 PS, Eq. 1.55 FS) in the jj dilute program

Reference: `qed3int_v2-14.pdf` Sec. 1 (Eqs. 1.16-1.55). Add the spacetime-averaged parity-symmetric and
flavor-symmetric condensates as cheap by-products of the existing diluted solve in
`jj_corr_dilute_claude.cu`.

## The operators
- PS (1.23): $\sigma_{PS}=\eta^\dagger\xi+\xi^\dagger\eta=\bar\Psi\Psi$.
- FS (1.55): $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta$ (massless $D_{ov}$ in the
  $(1-D_{ov}^\dagger)$ factor; $=\eta^\dagger\xi+\xi^\dagger V^\dagger\eta$, $V=D_{ov}-1$).

(Per flavor pair, $N_f=2$; spin index contracted.) **Spacetime average uses the simplicial site-area
measure $A_n=$ `dual_areas[n]`** (the sphere is irregular): $\langle\bar\sigma\rangle=\frac{\sum_{n,t}A_n\langle\sigma(n,t)\rangle}{N_t\sum_n A_n}$.
So the relevant trace is the AREA-WEIGHTED $\mathrm{tr}[A\,M]=\sum_{n,t}A_n M((n,t),(n,t))$ (a diagonal
$A=\mathrm{diag}(A_n)$ inserted before the trace), NOT the uniform $\mathrm{tr}\,M$. Stochastically
$\mathrm{tr}[A M]\approx (A\eta)^\dagger M\eta=\eta_A^\dagger\phi'$ with the area-weighted source
$\eta_A(t,n,i)=A_n\eta(t,n,i)$ ($A$ real -> $\eta_A^\dagger=\eta^\dagger A$). Denominator for the average:
$N_t\sum_n A_n$ ($=N_t\cdot4\pi$ on the unit sphere).

## Fermion contractions -> traces of $D_m^{-1}$
The single-pair action (Eq. 2.1, soft mass $m$) is $S=\eta^\dagger D_m\xi+\xi^\dagger D_m^\dagger\eta$,
$D_m=D_{ov}+m$. Writing it as $\big(\eta^\dagger\;\xi^\dagger\big)\!\begin{pmatrix}0&D_m\\ D_m^\dagger&0\end{pmatrix}\!\big(\eta\;\xi\big)^T$,
the inverse is $\begin{pmatrix}0&D_m^{-\dagger}\\ D_m^{-1}&0\end{pmatrix}$, so the only nonzero contractions are
$$
\langle\xi_a\,\eta^\dagger_b\rangle=(D_m^{-1})_{ab},\qquad
\langle\eta_a\,\xi^\dagger_b\rangle=(D_m^{-\dagger})_{ab},\qquad
\langle\xi\xi^\dagger\rangle=\langle\eta\eta^\dagger\rangle=0 .
$$
With the Grassmann reorder sign $\langle\eta^\dagger_b\xi_a\rangle=-(D_m^{-1})_{ab}$:
$$
\langle\eta^\dagger\xi\rangle=-\operatorname{tr}D_m^{-1},\qquad
\langle\xi^\dagger\eta\rangle=-\operatorname{tr}D_m^{-\dagger}=-\overline{\operatorname{tr}D_m^{-1}} .
$$
Hence
$$
\langle\sigma_{PS}\rangle=-\operatorname{tr}D_m^{-1}-\overline{\operatorname{tr}D_m^{-1}}
                        =-2\,\mathrm{Re}\,\operatorname{tr}D_m^{-1},
$$
$$
\langle\sigma_{FS}\rangle=-\operatorname{tr}D_m^{-1}
   +\operatorname{tr}\!\big[(1-D_{ov}^\dagger)D_m^{-\dagger}\big]
   =-\operatorname{tr}D_m^{-1}+\overline{\operatorname{tr}\!\big[(1-D_{ov})D_m^{-1}\big]} .
$$
(Sign of the $\xi^\dagger(1-D_{ov}^\dagger)\eta$ term: $\langle\xi^\dagger_a(1-D_{ov}^\dagger)_{ab}\eta_b\rangle
=-\operatorname{tr}[(1-D_{ov}^\dagger)D_m^{-\dagger}]$, and (1.55) subtracts it -> $+$ above.)

## Mass-pattern dependence of the propagator (Eqs. 3.50/3.51, 3.57/3.58, 3.60/3.61)
The forward contraction $\overline{\xi\,\eta^\dagger}$ is $D_m^{-1}$ in ALL cases (massless $D_{ov}^{-1}$ 3.50;
$m_F$: $D_{m_F}^{-1}$ 3.57; $m_P$: $D_{m_P}^{-1}$ 3.60). The BACKWARD contraction $\overline{\eta\,\xi^\dagger}$
differs:
- massless (3.51): $(D_{ov}^{-1})^\dagger$;  $m_F$ (3.58): $(D_{m_F}^{-1})^\dagger$  -> just the dagger, so
  $\langle\xi^\dagger\eta\rangle=\overline{\langle\eta^\dagger\xi\rangle}$ (conjugate; NO separate solve).
- $m_P$ (3.61): $(1+m_P)^{-1}(\tilde D_{m_P}^{-1})^\dagger$, $\tilde D_{m_P}=D_{ov}+m_P/(1-m_P)$ (3.62) --
  a DIFFERENT operator + the $(1+m_P)^{-1}$ factor; NOT a conjugate. In the dilute this $\tilde D$ propagator
  is already computed per pattern as `phimm` $=\tilde D_{m_P}^{-\dagger}\eta$ (parity (-) sink leg) -> reuse it,
  still no NEW solve.

=> THREE stored bilinears: `etadag_xi` $=\langle\eta^\dagger\xi\rangle$, `xidag_eta` $=\langle\xi^\dagger\eta\rangle$,
`xidag_1mDdag_eta` $=\langle\xi^\dagger(1-D_{ov}^\dagger)\eta\rangle$. Then
$\sigma_{PS}=$ `etadag_xi`$+$`xidag_eta` (1.23); $\sigma_{FS}=$ `etadag_xi`$-$`xidag_1mDdag_eta` (1.55).
Use the MASSIVE propagator with the SAME mass as the sea quarks (run `--mass-re/--mass-im` = sea mass; the
condensate reuses $D_m$/$\tilde D_{m_P}$ so it is automatically correct).

## Stochastic estimator -- reuse the existing $\phi'=D_m^{-1}\eta$ (and $\tilde D^{-\dagger}\eta$ for $m_P$)
The dilute already solves $\phi'=D_m^{-1}\eta$ once per dilution pattern. With $\langle\eta\eta^\dagger\rangle=I$:
$$
\operatorname{tr}D_m^{-1}\approx\eta^\dagger\phi',\qquad
\operatorname{tr}\!\big[(1-D_{ov})D_m^{-1}\big]\approx\eta^\dagger(1-D_{ov})\phi'
   =\eta^\dagger\phi'-\eta^\dagger(D_{ov}\phi') .
$$
So per pattern we need ONLY: the inner product $\eta^\dagger\phi'$ (already-available buffers) plus, for FS,
ONE extra $D_{ov}$ apply on $\phi'$ (no solve). PS and FS share $\eta^\dagger\phi'$. Cost = negligible vs a
solve; pure by-product of the jj solves. No K-apply, no current insertion, no $t_0$ origins (it is not a
source-origin object -- independent of the master-field superposition).

## Q1: reuse the solves? YES.
The condensates are traces of $D_m^{-1}$, obtained from the dilute's existing $\phi'$. Add a contraction
block right after the $\phi'$ solve (before/independent of the tp/sp/local/axial blocks). No extra solve.

## Q2: does time dilution help? YES (and it is FREE here).
The trace estimator's variance is $\sim\sum_{x\neq y}|D_m^{-1}(x,y)|^2$ (off-diagonal contamination; the
signal is the spacetime diagonal $\sum_x D_m^{-1}(x,x)$). Dilution into classes $\eta=\sum_d\eta_d$,
$\phi'_d=D_m^{-1}\eta_d$, replaces $\eta^\dagger\phi'$ by $\sum_d\eta_d^\dagger\phi'_d$, which DELETES every
off-diagonal pair $(x,y)$ whose source classes differ -> removes that part of the variance exactly.
- SPIN dilution (already ON): removes the $\sigma$ off-diagonal noise (2-spin trace).
- TIME dilution: removes off-diagonal pairs across time-classes. The dilute uses INTERLEAVED classes
  (`valence_claude.h time_spin_dilution`: interval $=N_t/t_{block}$, fill $t=t_s,t_s+$interval,$\dots$),
  so `td`$\ge2$ removes the dominant SMALL-$|\Delta t|$ time off-diagonals (e.g. td=2 kills all odd-$\Delta t$
  pairs incl. nearest-neighbor) -- exactly where $|D_m^{-1}|$ is largest. So time dilution is well-matched
  to the condensate.
- KEY: the dilute ALREADY runs the spin/time dilution patterns for jj and produces each $\phi'_d$. The
  condensate accumulated as $\sum_d\eta_d^\dagger\phi'_d$ over those SAME patterns inherits the
  variance reduction at ZERO extra solve cost. (Standalone, time dilution would cost $N_d\times$ solves;
  as a jj free-rider it is free.)

Caveat: dilution only reduces the STOCHASTIC variance; the gauge (config-to-config) fluctuation of this
disconnected single trace is separate and unaffected.

## Implementation chunk -- `jj_corr_dilute_claude.cu`  [DONE 2026-06-13]
Three stored BILINEARS (incl. Grassmann sign + per-case propagator), so the recipe gives (1.23)/(1.55):
- `etadag_xi`        $= \langle\eta^\dagger\xi\rangle = -\sum_d\eta_d^\dagger\phi'_d$  (ALL cases).
- `xidag_eta`        $= \langle\xi^\dagger\eta\rangle$: massless/$m_F$ $=\mathrm{conj}($`etadag_xi`$)$;
  $m_P$ $= -(1+m_P)^{-1}\sum_d\eta_d^\dagger\,$`phimm`$_d$.
- `xidag_1mDdag_eta` $= \langle\xi^\dagger(1-D_{ov}^\dagger)\eta\rangle$: massless/$m_F$
  $=-\overline{\sum_d\eta_d^\dagger(1-D_{ov})\phi'_d}$; $m_P$ $=-(1+m_P)^{-1}\sum_d\eta_d^\dagger(1-D_{ov}^\dagger)\,$`phimm`$_d$.

Code (in `jj_corr_dilute_claude.cu`):
1. Hit accumulators `acc_etadag_phi`, `acc_etadag_1mD_phi`, `acc_etadag_phimm`, `acc_etadag_1mDdag_phimm`
   (after `CsA`).
2. Per pattern after the $\phi'$ solve: `acc_etadag_phi += eta.dag(phi)` (all); if `!parity`,
   `acc_etadag_1mD_phi += eta.dag(op_oneMinusD(phi))`.  In the `if(parity)` block after the `phimm` solve:
   `acc_etadag_phimm += eta.dag(phimm)` and `acc_etadag_1mDdag_phimm += eta.dag(phimm) - eta.dag(op_DH(phimm))`
   (op_DH $=D_{ov}^\dagger$; $(1-D_{ov}^\dagger)$phimm).
3. Output: per-case `etadag_xi`/`xidag_eta`/`xidag_1mDdag_eta` -> `write_vec h0/condensate/{...}`.

Analysis (notebook): `sigma_PS = etadag_xi + xidag_eta`; `sigma_FS = etadag_xi - xidag_1mDdag_eta`;
spacetime-average by `/(n_sites*Nt)`; average over hits/configs.

Production: `run_jj_dilute_prod_claude.sh` time dilution turned ON (`TD=1`->`TD=2`) -> 4 patterns/hit, new
output dir suffix `_s1td2` (distinct from the in-progress `_s1td1`, no collision).

## Open questions
- OC1 (valence mass). Uses the simulation's $D_m$ (same soft mass as the dilute) for $\phi'$ -- chosen.
  Note FS's $(1-D_{ov})$ uses the MASSLESS overlap (`op_oneMinusD`/`M_D`) regardless.
- OC2 (normalization). Stored = AREA-WEIGHTED, full spacetime-summed traces $\mathrm{tr}[A\,M]$ (site
  weight $A_n=$ `dual_areas[n]` folded in via $\eta_A$; NOT divided). Spacetime AVERAGE is the analysis
  step: divide by $N_t\sum_n A_n$ ($=N_t\cdot4\pi$ unit sphere). Spin trace already included.
  Per-flavor-pair ($N_f=2$); apply any $N_f/2$ factor in analysis if wanted.
- OC3 (free check). DONE -- see "Free-field dense cross-check" below.

## Free-field dense cross-check  [program: `condensate_deter_free_claude.cu`]
The EXACT area-weighted traces are the stochastic estimator with a COMPLETE unit-source basis $\{e_x\}$
(the $N_{\rm hits}\to\infty$ limit; $\sum_x e_x e_x^\dagger=\mathbb{1}$ exactly), reusing the IDENTICAL
overlap operators (same `OverlapWMass`, `op_Dm`, `op_oneMinusD`, `op_DH`, `op_tilDm`), so it is
apples-to-apples with the dilute. For $x=0..N-1$ (free $U=1$):
- $\psi_x=D_m^{-1}e_x$ (op_DmH + op_Dmsq.solve); accumulate $A_{n(x)}\,\psi_x[x]\to\mathrm{tr}[A D_m^{-1}]$,
  and $A_{n(x)}\,[(1-D_{ov})\psi_x]_x\to\mathrm{tr}[A(1-D_{ov})D_m^{-1}]$ (massless/$m_F$ leg).
- $m_P$: $\phi^{mm}_x=\tilde D_{m_P}^{-\dagger}e_x$ (op_tilDm + op_tilDmsq.solve); accumulate
  $A_{n(x)}\,\phi^{mm}_x[x]$ and $A_{n(x)}\,[(1-D_{ov}^\dagger)\phi^{mm}_x]_x$.
- $n(x)=(x\bmod N_x)/N_S$, $N_x=N_S\,n_{\rm sites}$; $A_n=$ `dual_areas[n]`.
- SPEEDUP: free field is $t$-translation-invariant, so $D_m^{-1}(x,x)$ depends only on the spatial site ->
  sweep ONLY the $t=0$ block ($x\in[0,N_x)$, $N_x=24$ at L1) and multiply the traces by $N_t$ ($O(N_x)$
  solves, ~15-30 s vs ~30-50 min). A $t=1$ diagonal is checked against $t=0$ (must agree to ~CG tol) to
  confirm the invariance.
Same per-case combination as the dilute -> writes `condensate/{etadag_xi,xidag_eta,xidag_1mDdag_eta}` plus
`Vst` $=N_t\sum_n A_n$, to `data_<free>/condensate_deter_L<L>/cond.h5`. $N$ solves ($2N$ for $m_P$);
L=1 minutes. Runner `tmp_condensate_deter_claude.sh` (build to DISTINCT binary). Validation notebook
`condensate_validate_claude.ipynb`: jackknife the free stochastic dilute run (`h0/condensate/*` over hits),
spacetime-average both by `/Vst`, check stoch $\to$ dense for PS and FS (massless first; $m_F$, $m_P$ by
re-running both at that mass).
