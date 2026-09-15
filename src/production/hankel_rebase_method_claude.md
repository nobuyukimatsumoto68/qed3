# Block-Hankel + rebase for the $0^{++}$ $\sigma^2$/two-meson GEVP (interacting QED3)

Method note for the connected, vacuum-free variational analysis of the $0^{++}$ fermionic sector
($\sigma^2=(\bar\psi\psi)^2$ mixing with the two-meson state; the lattice precursor to the
Chester–Pufu $\sigma^2$–$F^2$ mixing, arXiv:1603.05582). Distillation on $S^2\times\mathbb R$, overlap
fermions. Driver: `gevp_twosigma_interacting_claude.py`.

## 1. Operators (interacting-safe)

- $\sigma^2$ — the target four-fermion operator $(\bar\psi\Phi\psi)^2$, $\Phi=$ Y00 ($\ell=0$) distillation smear.
- $O_A=\bar\psi\,\Phi\,\tilde\tau\,\psi$ — bilinear, $\tilde\tau=\tau(t,t)-\tfrac12$ (equal-time perambulator,
  GW contact $\tfrac12$ subtracted; `CONTACT=0.5`).
- $O_{2\sigma}(t)=\sigma_{00}(t)\,\sigma_{00}(t{+}\delta)$ — **time-split two-$\sigma$**, the genuine
  four-fermion two-meson interpolator ($\delta=$ `SPLIT`, default 1). A *bilinear* cannot reach the
  four-fermion two-meson (free-theory particle number), so the split is essential.
- **DROPPED for the interacting theory:** $O_{22}=\bar\psi Q_2\tilde\tau\psi$. The gauge field broadens
  the $S^2$ shells past their spacing (no 4/8/12 clustering at gsq 0.5/1.0/1.5), so the $\lambda=2$
  spectral projector $Q_2$ is ill-defined. (Fine in the free theory.)

## 2. Vacuum subtraction by diagram selection (NOT by fitting a constant)

$\langle\sigma^2\sigma^2\rangle$ has 10 Wick diagrams (`diags_pair`, labels A..J). The **connected** ones
**A,B,C,D,E,G** decay; **F,H,I,J** are disconnected tadpole/self-loop products that sit at a (slowly
varying) constant — the condensate/vacuum piece. Key findings (`diag_corr_linear_claude.py`,
`diag_BE_vs_2mPS_claude.py`, `diag_sum_vs_mPS_claude.py`):

- **B+E** $=C_S^2-T_S$ (direct + crossed) $\to$ the two-meson, effmass $\approx2m_{PS}+O(0.05)$ (weak scattering shift).
- **A,C,D,G** carry a one-meson-rich tower (effmass $\approx0.46$, between $m_{PS}$ and $2m_{PS}$).
- **F,H,I,J** = flat vacuum constants; a "double plateau subtraction" of the full sum does NOT clean up
  (F,H,I leave a slow $\sim0.1$ residual + variance noise). **Exclusion beats subtraction.**

So the vacuum-free correlator is **$\langle\sigma^2\sigma^2\rangle_{\rm conn}$ = A+B+C+D+E+G** (CONN_IDX).
$\langle O_A\sigma^2\rangle$ reduces to a single connected **triangle** ($=C_{2A}$; Semi/Disc $\propto\langle\sigma\rangle\approx0$
after contact subtraction; OA_Dp is vacuum). $\langle O_{2\sigma}\cdot\rangle$ built from the single-$\sigma$
propagator $M(t,s)=-\mathrm{Tr}[\Phi(t)\tau(t,s)\Phi(s)\tau(s,t)]$.

**Connected GEVP matrix, NO identity** (vacuum-free by construction; `CONN=1`, `NOOA=1` for the 2-op):
$$
C(t)=\begin{pmatrix}\text{A..G}(t) & C_{2s}(t)\\ C_{2s}(t) & C_{ss}(t)\end{pmatrix}\ \{\sigma^2,O_{2\sigma}\},
\quad\text{or the }3\times3\ \{\sigma^2,O_A,O_{2\sigma}\}.
$$
Energies are invariant under per-operator rescaling, so the (factor-2) normalization convention is irrelevant.

## 3. Block-Hankel proliferation + rebase (the trick)

Enrich the basis with time-shifts (generalized pencil of functions; Aubin–Orginos 1010.0202) WITHOUT new
operators, then variationally truncate ("rebase") onto the leading states:

1. **Block-Hankel:** $\hat C(t)_{(a,i),(b,j)}=C\big(t+(a{+}b)\,\text{SHIFT}\big)$, $a,b=0..\text{NSH}-1$.
   Naming: "Dt=1,2" $\Rightarrow$ NSH=3, SHIFT=1 (blocks at offsets 0,1,2); "Dt=3 only" $\Rightarrow$
   NSH=2, SHIFT=3 (0,3); "Dt=4,8" $\Rightarrow$ NSH=3, SHIFT=4 (0,4,8). Costs range: $t_\max=$ twin$-2(\text{NSH}-1)\text{SHIFT}$.
2. **Rebase:** solve the enlarged GEVP at $t_{\rm reb}$ (metric $t_0^{\rm reb}$), keep the leading NKEEP
   eigenvectors $V$ (largest $\lambda$ = lightest, needs $t_{\rm reb}>t_0^{\rm reb}$), project
   $C_{\rm reb}(t)=V^\top\hat C(t)\,V$ (NKEEP$\times$NKEEP), and read the effmass from a final GEVP on $C_{\rm reb}$.
3. **$Z$ overlaps** (Blossier 0902.1265): $Z_n^a=(C(t_0)v_n)_a\,e^{E_n t_0/2}$, $v_n$ normalized
   $v_n^\top C(t_0)v_n=1$ (`zmatrix`, heatmaps).

**Why it works:** the time-shifts fold in the *well-overlapped, low-contamination large-$t$* information
while the small-$t$ data supplies statistics; rebasing then locks in the clean $N$-state subspace. Net:
plateaus reached much earlier (e.g. $m_1$ at $t=5$ drops $0.96\to0.67$) at balanced S/N.

## 4. Knobs, findings, and warnings

- **`SHIFT` / `NSH` (which Dt's):** larger shifts kill contamination faster but eat range and add noise to
  the light state. **Dt=1,2 (or Dt=3) is the sweet spot**; wider spacings ($[0,3,6]$, $[0,4,8]$) reach a
  *lower* clean $m_1$ plateau (see the variational note below — lower is *better*, not biased).
- **BLOCK-COUNT CEILING (2-op base = 3 blocks, NSH=3).** The connected 2-op signal is rank-2, so a Hankel of
  $n_b$ blocks is $2n_b\times2n_b$ but carries rank $\sim2$. Empirically **3 blocks ($6\times6$) is the maximum
  that stays conditioned** ($[0,1,2]$, $[0,2,4]$, $[0,3,6]$, $[0,4,8]$ all fine); **4+ blocks are
  rank-degenerate** ($[0,2,4,6]$ $8\times8$ and $[0,1,\dots,6]$ $14\times14$ give NaNs / huge errors — the
  $t_0$ metric goes near-singular), regardless of shift spacing. To go wider you need MORE base operators, not
  more time-offsets.
- **VARIATIONAL DIRECTION (crucial):** GEVP effmasses approach the true energy **from above** (upper bound),
  so a cleanly-reached *lower* plateau is the *better* estimate — NOT a downward "over-cleaning bias". A
  narrow-shift value read on a still-descending shoulder overestimates; the wide-shift flat plateau is closer
  to truth. (Earlier "low-bias artifact" framing was wrong and is retracted.)
- **`REBT` (rebase time), `REBT0` (rebase metric):** for a **2-op base the rebase point is irrelevant** —
  the $4\times4$ Hankel is rank-2, so the leading-2 subspace is fixed ($t_{\rm reb}=1..9$ gave *identical*
  $m_0=0.462$, $m_1=0.646$). Rebase choice only matters once the base is 3+ ops (NKEEP < available states).
  $t_{\rm reb}=0$ is the degenerate metric point ($t_{\rm reb}\le t_0^{\rm reb}$ flips light/heavy).
- **`NKEEP` (rebased states):** truncating to fewer states than the basis supports is where the real
  variational gain (and risk) lives.
- **WARNING — tilted state-1:** some knob choices leave the excited/two-meson effmass *sloped* rather than
  flat — a bad signal (residual contamination or a compromised subspace). A flat $m_1$ is the acceptance
  criterion. "Dt=1,2 rebase 2@$t{=}5$" gives the best balance (fast plateau, flat $m_1$, good S/N).
- **Autocorrelation:** naive (binsize-1) jackknife underestimates errors; a binsize scan
  (`BINSCAN=1`) shows the error plateaus by block $\approx8$–16 (ground $\sim1.7$–$2\times$, two-meson
  $\sim1.1$–$1.4\times$). Quote errors at **`BINSIZE`$\approx$10–16**.
- **Multiple rebasing** (rebase at several $t$'s / staged truncation) is unexplored — a natural extension.

## 5. Result (Nf2 gsq1.0 L1, 400 cfg, binsize 10)

Preferred config: **offsets $[0,3,6]$, 2-op, $T_0=3$, rebase 2@$t{=}4$** (clean flat plateaus). Correlated
constant fits (jackknife covariance; $n_{\rm bin}=40$):
$$
m_0=0.4296(166)\ (t\in[4,6],\ \chi^2/\text{dof}=0.06/2),\qquad
m_1=0.6145(58)\ (t\in[12,14],\ \chi^2/\text{dof}=0.07/2).
$$
(The narrow-shift Dt=1,2 pick reads $m_1\approx0.644$, but on a still-descending shoulder — the $[0,3,6]$
plateau is the clean, lower, better value; see the variational note.)

### 5a. SUB-THRESHOLD (candidate BOUND) two-meson — the equal-footing test
Extract $m_{PS}$ with the SAME Hankel machinery so the threshold $2m_{PS}$ is on equal footing with $m_1$:
1-op base $M(t)=-\mathrm{Tr}[\Phi\tau\Phi\tau]$, offsets $[0,3,6]$ (GPOF resolves a small tower), $T_0=3$,
rebase 2@$t{=}4$, state 0 $=m_{PS}$. Ground plateau $m_{PS}\approx0.32$, so $2m_{PS}(\text{clean})=0.6387(45)$.
Then $m_1=0.6145(58)$ sits **below $2m_{PS}$ in every $m_{PS}$ window** ($\Delta=m_1-2m_{PS}\approx-0.014$ to
$-0.027$, $1.4$–$3.3\sigma$; sign robust). This is the trustworthy binding evidence — a clean-plateau $m_1$
compared to a shoulder $m_{PS}$ would be apples-to-oranges. **CAVEAT:** $m_{PS}$ itself drifts slowly
($0.322\to0.314$), so $|\Delta|$ carries a window systematic, but even the most conservative $2m_{PS}=0.628$
keeps $m_1$ below. Consistency with the literature is unresolved (this is $N_f$-dependent; the project's own
no-SSB / near-conformal finding at $N_f{=}2$ makes a near-threshold two-meson the natural expectation and deep
binding a surprise). To call it a genuine bound state rather than a finite-volume (Lüscher) level shift needs a
finite-volume threshold analysis and repetition at other $g^2$/$N_f$/$L$. Files: `hankel_rebase_dt036_fit_claude.py`
(the fit), `hankel_mPS_dt036_claude.py` (equal-footing $2m_{PS}$), `hankel_rebase_nohankel_plot_claude.py`
(no-Hankel baseline: both states plateau only at large $t$).

## 6. Files

- `gevp_twosigma_interacting_claude.py` — driver. Env knobs: `CONN NOOA SPLIT T0 BINSIZE BINSCAN
  HANKEL NSH SHIFT REBT NKEEP REBSCAN RREAD REBT0 NCFG`. Produces effmass, $Z$ heatmap, binscan, rebase-scan figs.
- Diagram tools: `diag_effmass_claude.py` (10-diagram effmass, `diags_pair`), `diag_corr_linear_claude.py`
  (per-diagram linear corr, panels), `diag_BE_vs_2mPS_claude.py`, `diag_sum_vs_mPS_claude.py` (`DIAGS=`),
  `diag_OA_sig2_linear_claude.py` (⟨$O_A\sigma^2$⟩ diagrams).
- Free-theory validation (where $O_{22}$/$Q_2$ works): `gevp_twosigma_split_free_claude.py`,
  `gevp_timesplit_free_claude.py` (GPOF), `oa_o22_shell_decomp_free_claude.py`.

Refs: distillation Peardon et al. 0905.2160; GPOF Aubin–Orginos 1010.0202; GEVP $Z$-factors Blossier
et al. 0902.1265; mixing Chester–Pufu 1603.05582.
