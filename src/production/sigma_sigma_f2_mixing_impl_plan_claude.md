# $\sigma^2$ / vector$^2$ -- $F^2$ mixing in the $0^{++}$ channel: implementation plan

## Physics goal

Measure, on the QED3 lattice, the large-$N_f$ CFT prediction that the $0^{++}$ glueball $F^2$
mixes with the two-scalar (four-fermion) singlet $(\bar\psi\psi)^2$. Build a coupled $0^{++}$
GEVP in the operator basis $\{F^2\ \text{shapes},\ \text{two-meson }\sigma\sigma,\ \text{vector}^2\}$,
extract the eigenstates, and check that the light branch is mostly $F^2$ with $\Delta_-\to 4$
(not 3) as $N_f$ grows.

**Sources (mandatory citation):**
- S. M. Chester, S. S. Pufu, arXiv:1603.05582 -- the $0^{++}$ two-operator mixing
  $\mathcal{O}_1=(\bar\psi\psi)^2$, $\mathcal{O}_2=F^2$ (Eqs. 70-93); $\Delta_-=3.30,3.65,3.77$ at
  $N=N_f=2,4,6$ (two-component; this code is two-component so $N=N_f$).
- In-repo `../../t_orthog_largeNf_claude.md` (the large-$N_f$ thread). Also records that the
  singlet vector$^2$ $\mathcal{O}_3=(\bar\psi\gamma_\mu\psi)^2$ is the redundant/gauged (EOM)
  operator, $\langle\mathcal{O}_3\mathcal{O}_3\rangle=0$ at LO only.

## Operator content (corrected picture, NM 2026-08-23)

### Scalar densities and parity
$\sigma_{PS},\sigma_{FS}$ are FLAVOR-SINGLET, flavor-DIAGONAL bilinears (the $\xi,\eta$ are the
within-flavor GW/chiral components of one Dirac fermion, flavor index summed) -- so their
zero-momentum loops are NONZERO. (Earlier "scalar disc = 0 because off-diagonal" was WRONG.)
- $\sigma_{PS}=\eta^\dagger\xi+\xi^\dagger\eta$ -- **parity-even**, definite parity.
- $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta$ -- **no definite parity on the
  lattice** (the $(1-D_{ov}^\dagger)$ GW insertion mixes parities; sharpens only as $a\to0$).

Zero-momentum ($l=0$) loops, per timeslice, from one inversion $\phi=D_m^{-1}\eta$:
$$
J_1(t)=\mathrm{tr}[W_0 D_m^{-1}]_t,\qquad J_{1mD}(t)=\mathrm{tr}[W_0(1-D_{ov})D_m^{-1}]_t,\qquad W_0=l{=}0\ \text{vertex},
$$
$$
J_{PS}(t)=2\,\mathrm{Re}\,J_1(t),\qquad J_{FS}(t)=J_1(t)-\overline{J_{1mD}(t)}.
$$
Both share the single solve $\phi$ ($J_{1mD}$'s $(1-D_{ov})\phi$ is a mat-vec, no extra solve).

### Two-meson operators -- FULL pairing basis (no parity pruning)
$$
\mathcal{O}_{ab}(t)=J_a(t)J_b(t),\qquad a,b\in\{FS,PS\}:\quad \{FS{\cdot}FS,\ PS{\cdot}PS,\ FS{\cdot}PS\}.
$$
Because FS lacks definite parity at finite $a$, EVERY pairing carries a $0^{++}$ component --
there is no parity veto to drop any of them (it re-emerges only as $a\to0$). Parity is still a
symmetry of the ACTION, so the correlators self-project: $\langle F^2\,\mathcal{O}_{ab}\rangle$
picks out the parity-even ($0^{++}$) content automatically; the parity-odd part drops. CP's
$(\bar\psi\psi)^2$ ($0^{++}$) is the GW/continuum combination inside this basis (mostly the
FS-containing pairings, since PS is the clean parity-even one); the GEVP isolates it.

### Vector$^2$ -- redundant channel, FREE from existing 2-hit data
$\mathcal{O}_3(t)=\sum_x J_\mu J^\mu(x,t)=\sum_{a,l,m}J_{a,l,m}(t)\overline{J_{a,l,m}(t)}$ from the
existing singlet-vector disc loops `h0/disc/ylm/s{1,2,3}/l/m/J` (dir `corr_ylm_disc_tb2`, 2 hits).
By EOM $J^\text{singlet}_\mu\propto dF$ is redundant, so $\langle\mathcal{O}_3\mathcal{O}_3\rangle$
vanishes at large-$N_f$ LO (CP Fig. 11) but is $O(1/N_f)$-nonzero at our $N_f$. Near $F^2$ in the
GEVP (it is $\propto (dF)^2$); the finite-$N_f$ splitting is the signal. Signal check underway
(`vec2_signal_claude.py`).

### $F^2$ (gauge)
$l{=}0$ $F^2$ shape operators from `glue_f2_v2_shapes` (14 shapes: 7 $F^2$ + 7 $F^4$).

## Reduction: all correlators are products of per-timeslice loop series

$C_{xy}(t)=\langle J_x(t)J_y(0)\rangle_c$, $x,y\in\{FS,PS\}$ (single-loop 2pt). Then:
$$
\langle \mathcal{O}_{ab}(t)\mathcal{O}_{cd}(0)\rangle_{\text{disc}}=C_{ac}(t)C_{bd}(t)+C_{ad}(t)C_{bc}(t),
\qquad
\langle F^2(t)\,\mathcal{O}_{cd}(0)\rangle=\langle O_F(t)J_c(0)J_d(0)\rangle_c .
$$
So the disconnected two-meson block, the $F^2$ mixing, and the vector$^2$ block are ALL built
by forming per-config products of the elementary series $\{J_{FS}(t),J_{PS}(t),J_{a,l,m}(t),O_F(t)\}$
and gauge-averaging. No four-point inversion for any of this. (The genuinely CONNECTED four-point
of $\langle\mathcal{O}_{ab}\mathcal{O}_{cd}\rangle$ is the deferred valence piece -- next discussion.)

## Measurement: combined $l{=}0$ scalar-loop driver (stride-10, nhits=7, mrhs)

New driver `jj_sigma_loops_stoch_claude.cu` in `src/production/` (port + trim of
`src/both_3d/jj_local_ylm_scalar_disc_stoch_claude.cu`): per config, per hit, $l{=}0$ only,
compute $J_1(t)$ and $J_{1mD}(t)$ (=> $J_{FS},J_{PS}$), sharing the inversion.

### Statistics / cost (why stride-10, nhits=7)
Two-meson operators are products of two loops -> unbiased only from DISTINCT hits ($\binom{H}{2}$ pairs).

| scheme | configs | hits | inversions | pairs/cfg | total pairs |
|--------|---------|------|-----------|-----------|-------------|
| stride-2, nhits=2  | $5N$ | 2 | $10N$ | 1  | $5N$  |
| stride-10, nhits=7 | $N$  | 7 | $7N$  | 21 | $21N$ |

stride-10/nhits=7 is cheaper AND gives more unbiased estimates. Caveat: the 21 pairs cut only the
STOCHASTIC variance ($\sim 2/[H(H-1)]$ per config); the GAUGE variance is set by config count, and
stride-10 samples $5\times$ sparser (partly offset by less autocorrelation). Net win iff the loop
is stochastic-noise-dominated -- validate empirically (the vec2 / gate checks indicate the balance).

### Unbiased estimator
Store all per-hit loops; pair distinct hits in analysis.
$\mathcal{O}_{ab}(t)=\frac{1}{H(H-1)}\sum_{h\ne h'}J_a^{(h)}(t)J_b^{(h')}(t)$;
$\langle F^2\mathcal{O}_{cd}\rangle$ uses distinct hits for the two loop factors ($O_F$ exact).
tb2 time-dilution handles the same-hit cross-time contamination for separated $t$.

### mrhs blocking (solver only; h5 format UNCHANGED)
Block the 7 stochastic RHS into one `BlockedMat` solve (the C6f / GRAD_L4 multi-RHS machinery),
extract the 7 solutions, form 7 per-hit $l{=}0$ loops, write them as the SAME per-hit datasets.
Only the CG call changes; expect ~2-3x on the solve on top of the cheaper inversion count.

### Driver / h5
- keys `h0/sigma_loops/s0/J`=$J_1(t)$, `h0/sigma_loops/s0_1mD/J`=$J_{1mD}(t)$, complex len $N_t$, PER HIT.
- $a_t$ auto-derived from ens-dir (reuse `at_from_ensdir`).
- resumable / `complete`-gated per (config,hit); `--nhits` default 7; NO rm / NO kill in scripts.
- NEW output subdir (e.g. `corr_sigma_loops_nhits7`), distinct from the vector `corr_ylm_disc_tb2`.

## $F^2$ per-timeslice dump (gauge side)
Add ONE dataset `O` (nops$\times N_t$, the in-memory `obs`) to `glue_f2_v2_shapes_claude.cu` at the
h5 write block (~line 493); config-match to the loop driver. Keep `F_corr_blk`/`F` for the standalone
$F^2$ GEVP. Gauge-only, cheap.

## Analysis: coupled $0^{++}$ GEVP
`sigma_f2_gevp_analysis_claude.cu` (mirror `glue_gevp_analysis_claude.cu`, streaming binned jackknife):
build the correlation matrix over $\{F^2\ \text{shapes},\ FS{\cdot}FS,\ PS{\cdot}PS,\ FS{\cdot}PS,\ \text{vector}^2\}$
from unbiased per-config products, vacuum-subtract, $t_0$-metric GEVP,
$\Delta_\text{eff}=-\log\lambda/(dt\,a_t)$, eigenvector content ($F^2$ vs four-fermion). Compare
$\Delta_-$ to CP $\{3.30,3.65,3.77\}$.

## GATING CHECK (do first): $m{=}0$ disconnected susceptibility
Before the full driver, confirm $\langle J_{FS}(t)J_{FS}(0)\rangle_c$ and $\langle J_{PS}\cdots\rangle_c$
are NONZERO at exact $m{=}0$ (one-point $\langle J\rangle=0$ by no-SSB, but the connected loop-loop
piece is generically alive -- the $\eta'$/topological-susceptibility-type object). Check no exact-chiral
Ward identity flattens it. If it vanishes, the scalar $\sigma^2$ channel is dead and only vector$^2$
survives. (vec2 signal check is the vector analog, running now.)

## Connected four-point: solve strategy (Ch.5 of `../../qed3int_v3-4.pdf`; supersedes v3-2)

**v3-4 corrections (NM 2026-08-28) -- authoritative:**
- $V_S$ (5.18): $(\Pi_0,\Pi_t)$, 3 solves, ONE $\Pi_t$ -> sliceable single one-end trick. Diagrams C, D.
- $S_S$ (5.19): $(\Pi_0,\Pi_t,\Pi_t)$, **4 solves, TWO $\Pi_t$**. $\mathfrak{A}_A=-S_S$.
- $T_S$ (5.20): $(\Pi_t,\Pi_0,\Pi_t)$, **4 solves, TWO $\Pi_t$**. $\mathfrak{A}_B=-T_S$ (NOT 5 solves, NOT unused --
  my earlier v3-2 read was wrong on both). A and B BOTH needed.
- TWO $\Pi_t$ = an internal "$t\to t$ return" a single slice can't absorb = the **all-to-all propagator problem**.
  **INTENDED METHOD (NM 2026-08-29): DISTILLATION / LapH** -- perambulator $\tau_{t',t}=V^\dagger D_{ov}^{-1}V$
  ($V$ = low $U(1)$-covariant $S^2$-Laplacian eigenvectors) gives the all-to-all once; internal $\Pi_t\to\tau(t,t)$,
  every diagram = a perambulator trace, AND it de-noises the disc loops. Full writeup:
  **`distillation_for_two_meson_claude.md`**. (Alt/lightweight = NM's two-one-end-tricks, source at 0 AND t,
  bespoke per topology -- superseded by distillation as the general fix.) Conn driver DEFERRED (mixing CHECK
  first). $C_S$/E/F/G unaffected (single free $\Pi_t$).


The genuinely connected four-point of the SAME scalar (the two-meson operator's own conn correlator,
needed for the $\{FS{\cdot}FS,PS{\cdot}PS,FS{\cdot}PS\}$ diagonal blocks of the GEVP; the $F^2$ mixing
itself is disc). Ten topologies A..J (Fig.1), weights $\{4,2,4,4,2,1,4,1,1,1\}$;
$G_4=\langle S_4\rangle_g+\langle\tilde S_4\rangle_g$ (the $S$ and $\tilde S$ vertices never mix in a
trace). Massless propagator $G=D_{ov}^{-1}$. FS vertex $\tilde S=-(1-D_{ov}^\dagger)$ is a $D_{ov}^\dagger$
APPLY (mat-vec), not a backward inverse; in overlap each apply is one sign-function eval (inner solve).

### Six factors (all off one source $\phi_0$; $\Pi_\tau$ = timeslice-$\tau$ projector)
Forward solves off $\phi_0$:
$$
g_1=D_{ov}^{-1}\phi_0,\qquad g_2=D_{ov}^{-1}\,\Pi_t\,S\,g_1 .
$$
- $D_S(\tau)=\langle\langle\phi_\tau^\dagger S D_{ov}^{-1}\phi_\tau\rangle\rangle$ (5.15) -- plain $\ell{=}0$
  loop (= the disc-block loop $J$); **1 solve**.
- $D'_S(\tau)=\langle\langle\tilde\phi_\tau^\dagger S D_{ov}^{-1}\Pi_\tau S D_{ov}^{-1}\tilde\phi_\tau\rangle\rangle$
  (5.16) -- extended loop; **2 solves**.
- $C_S(0{\to}t)=\langle\langle\phi_0^\dagger S g_2\rangle\rangle$ (5.17) -- meson corr; **2 solves** ($g_1,g_2$).
- $V_S,S_S,T_S$ (5.18-5.20) -- nested chains; naively 3/4/5 solves, collapsed below.

### Meet-in-the-middle for $C_S,V_S,S_S,T_S$ (the reuse)
Fix the RIGHT half at the reusable second solve $g_2$; build the LEFT half by backward $D_{ov}^{-\dagger}$
solves off $\phi_0$ (genuine backward CG -- NO $\gamma_5$-herm, so $D_{ov}^{-\dagger}\ne\gamma_5 D_{ov}^{-1}\gamma_5$):
$$
L_V=D_{ov}^{-\dagger}S^\dagger\phi_0,\quad
L_S=D_{ov}^{-\dagger}S^\dagger\Pi_0 L_V,\quad
L_T=D_{ov}^{-\dagger}S^\dagger\Pi_t L_S\qquad(L_V\subset L_S\subset L_T,\ \text{nested}).
$$
Junctions (inner products, meeting point chosen so the right half reproduces $g_2$):
$$
V_S=L_V^\dagger\,\Pi_0\,S\,g_2\ (1{+}2),\qquad
S_S=L_S^\dagger\,\Pi_t\,S\,g_2\ (2{+}2),\qquad
T_S=L_T^\dagger\,\Pi_t\,S\,g_2\ (3{+}2).
$$
V is unbalanced (1 back, 2 fwd) because pinning the right half to $g_2$ (not balance) sets the split.
**Solve count: 2 forward ($g_1,g_2$) + 3 backward ($L_V,L_S,L_T$)** per source cover ALL of
$\{C_S,V_S,S_S,T_S\}$, vs $2{+}3{+}4{+}5=14$ done independently. $g_1$ is $S$-independent (shared PS/FS);
$g_2$ and the $L$'s are $S$-dependent (recompute per PS/FS, but shared across the four chains).

### $D_S$ + $D'_S$ combined loop driver -- 2 hits, time+spin diluted (NEW, massless)
$D_S$ has a MASSIVE implementation already; we need the MASSLESS one, measured TOGETHER with $D'_S$ so the
solves are shared. Per hit, ONE time+spin-diluted source $\phi$: first solve $u=D_{ov}^{-1}\phi$ gives
$D_S(\tau)=\langle\phi^\dagger S u\rangle_\tau$ directly AND is $D'_S$'s first solve; $D'_S$ adds a second
solve $D_{ov}^{-1}\Pi_\tau S u$ -> $D'_S(\tau)=\langle\phi^\dagger S(\cdot)\rangle_\tau$. So 2 solves/hit
yield both. **2 hits** (needed anyway for the $D_S$ squares H/I/J). With 2 hits, $\mathfrak{A}_F=D'_S(0)D'_S(t)$
uses DISTINCT hits -> fully unbiased for all $t$ (the earlier single-hit interval-$N_t/2$ bias at $t{=}N_t/2$
is moot; interval-$N_t/2$ stays only as dilution granularity).

RESOLVED (middle $\Pi_\tau$): NO per-timeslice second solve. Project the middle to the whole dilution CLASS
$P_\text{class}=\sum_{\tau\in\text{class}}\Pi_\tau$ and do ONE second solve per (class,spin): $\chi=D_{ov}^{-1}(P_\text{class}S\phi)$.
The $\tau'{=}\tau$ diagonal is $D'_S(\tau)$; $\tau'{\ne}\tau$ terms are propagator-suppressed by $\ge$ interval
(same suppression the $D_S$ one-point loop relies on) -> with interval $N_t/2$ the residual sits at $t{=}N_t/2$.

FS VERTEX (RESOLVED, NM 2026-08-28): $\tilde S=-(1-D_{ov}^\dagger)$ WITH dagger (Eq. 5.3). $D_S$ keeps the
one-point convention (store $J_1$ with $(1-D_{ov})$, FS $=J_1-\overline{J_{1mD}}$ at analysis). But $D'_S$ is a
TWO-vertex factor and conj does NOT recover it (conj flips $D^{-1}\!\to\!D^{-\dagger}$), so the FS extended loop
applies $(1-D_{ov}^\dagger)$ LITERALLY at both vertices (sign $(-1)^2{=}+$): middle $(1-D_{ov}^\dagger)\phi$,
project, solve, sink $(1-D_{ov}^\dagger)$. -> $D'_{1mD}$ is the literal FS factor, one piece, no mixed term.
Driver: `jj_sigma_loops_stoch_claude.cu` step (E) uses `op_oneMinusDdag`; step (C)/$D_S$ keeps `op_oneMinusD`.

### Square-bias handling (E,H,I,J only)
The $\langle\langle\cdot\rangle\rangle^2$ in (5.10/5.13/5.14) is the $H{\to}\infty$ object; at finite hits
the plain square $\bar X^2$ is biased by $+v/H$ (self-variance), NOT removed by more configs. FIX reuses
the SAME hits (no extra solve): distinct-hit product $\frac{1}{H(H-1)}\sum_{h\ne h'}X^{(h)}X^{(h')}$
($=X^{(1)}X^{(2)}$ at $H{=}2$). Squared objects: $C_S$ (in E), $D_S$ (in H,I,J).
- $D_S$: a NEW massless **2-hit** loop (7-hit DEFERRED), measured together with $D'_S$ (shared first solve);
  the square's two factors are the 2 distinct hits ($X^{(1)}X^{(2)}$) -> H/I/J unbiased (1 pair/config).
- $C_S$ (E): DO NOT need 2 hits in the new driver. $C_S(0{\to}t)$ IS the existing single-meson conn
  (`corr_ylm_conn_t00_nhits1_s1`, nhits=1, independent source). Pair it with the new driver's $C_S$
  (the $g_2$ right-half, an independent source) -> $E=C_S^{\text{existing}}\cdot C_S^{\text{new}}$ is
  unbiased (two independent single-hit estimates, NM 2026-08-28). Requires: SAME massless $C_S$
  object/normalization + config-match + independent RNG. So ALL 10 diagrams are unbiased at nhits=1.
- $D'_S,V_S,S_S,T_S$ single-hit-unbiased. Full derivation: `two_meson_square_bias_note_claude.md`.

### 10-diagram assembly (factors | #indep sources | bias-at-nhits1)
A(4)=$S_S$ | 1 | ok. B(2)=$S_S$-var | 1 | ok. C(4)=$D_S(t)\,V_S$ | 2 | ok. D(4)=t-rev C | 2 | ok.
E(2)=$C_S^2$ | 1 | ok via $C_S^{\text{existing}}\!\times\!C_S^{\text{new}}$ (INDEPENDENT source vectors -> distinct seed). F(1)=$D'_S(0)D'_S(t)$ | 2 (distinct hits -> unbiased all $t$) | ok.
G(4)=$D_S(0)D_S(t)C_S$ | 3 | ok. H(1)=$D_S(0)^2 D'_S(t)$ | 2 | ok ($D_S$ 7-hit). I=t-rev H | ok ($D_S$ 7-hit).
J(1)=$D_S(0)^2 D_S(t)^2$ | 2 | ok ($D_S$ 2-hit). ($D_S$ from the existing 2-hit disc loop [7-hit deferred];
E from existing$\times$new $C_S$ -> ALL 10 unbiased at nhits=1. KEY: new-driver $\phi_0$ seeded DIFFERENTLY
from the FNAL conn source AND the disc-loop source.)

## SCOPE CHANGE (NM 2026-08-29): mixing CHECK first, NOT the full GEVP
Goal downgraded to "do the scalars mix with $0^{++}$?" -- the cross-correlator
$$\langle F^2(t)\,\sigma\sigma(0)\rangle_c=\big\langle O_F(t)\,[\,D_S(0)^2+D'_S(0)\,]\big\rangle_c .$$
$\sigma\sigma(0)$ is a four-fermion op ALL on timeslice 0; its self-contraction there = $D_S(0)^2$ (disc, distinct
hits) + $D'_S(0)$ (connected single loop). Both correlate with $O_F(t)$ ONLY through the gauge field (no fermion
line 0->t), so this is **all-to-all-free** -- needs just the loop driver ($D_S$, $D'_S$), NOT the conn four-point
($S_S/T_S$ deferred). The $\sigma\sigma(t)\sigma\sigma(0)$ self-correlator (which needs the conn) is DEFERRED.
- LAUNCHED (2026-08-29): loop driver on L1 massless via `run_sigma_L1_claude.sh` (GPU1 MPS-2, stride-10, nhits=2).
- GLUE $O_F(t)$ GAP: glue driver writes only `F_corr_blk` (F^2 auto-corr) + `F` ($\sum_t O_F$), NOT the per-timeslice
  $O_F(t)$. For the $t$-resolved cross we must ADD an `O` dataset (nops x Nt) to `glue_f2_v2_shapes_claude.cu`
  (~L493, write `obs`) and re-run glue on the L1 configs (cheap). DEFERRED (NM tired of glue). Crude integrated
  susceptibility $\langle(\sum_t O_F)(\sum_t S)\rangle_c$ uses the EXISTING `F` -- fallback if we want a number now.

## Measurement plan for the conn/extended part (NM, 2026-08-28)
- **L=1 only** to begin (full analysis; higher L unlikely to show signal). stride-10.
- **Loop driver** (NEW, massless): $D_S$ + $D'_S$ combined, time+spin diluted, **2 hits**, shared first
  solve (2 solves/hit). 7-hit idea DEFERRED.
- **Conn driver**: $\{C_S,V_S,S_S,T_S\}$ meet-in-the-middle (2 fwd + 3 back), **nhits=1**; E pairs the new
  $C_S$ with the existing FNAL conn $C_S$.
  Both L=1, stride-10, **GPU1 + MPS 2-pack**. Measure COST + look at data, then decide on more statistics.
- **SEED INDEPENDENCE (critical).** The new conn/$D'_S$ sources must be independent of BOTH the disc loop
  set AND the existing FNAL conn, else products bias:
  (a) E ($C_S^2$) = new $C_S$ $\times$ existing FNAL conn $C_S$ -> new $\phi_0$ must differ from the FNAL conn seed;
  (b) C,D,G,H,I multiply a disc $D_S$ (7-hit set) by a new-conn factor ($V_S/C_S/D'_S$) -> new sources must
  differ from the disc-loop seed.
  Implement as DISTINCT seed namespaces/tags per stream: {disc-loop, conn-CVST, $D'_S$, FNAL-conn} all
  different. (J $=D_S^2 D_S^2$ is pure disc-loop, no conn factor -> assembled from the 7-hit set alone.)
- Drivers (to write): (1) loop `jj_sigma_loops_stoch_claude.cu` ($D_S$+$D'_S$, massless, 2 hits, shared
  solves) -- port the both_3d massive extended-loop; (2) conn `jj_sigma_conn_stoch_claude.cu` (C/V/S/T
  meet-in-the-middle, nhits=1); (3) handoff `run_sigma_L1_claude.sh` (user runs; GPU1 MPS-2; no rm/no kill;
  complete-gated; distinct seed namespaces per stream).

## Ordered chunks
1. **Gating check** -- $m{=}0$ disc susceptibility of $J_{FS},J_{PS}$ (small run or from a quick
   loop measurement on a few configs). Files: gate script/notebook.
2. **vector$^2$ signal check** -- `vec2_signal_claude.py` (existing data). RUNNING.
3. **Combined loop driver** -- `jj_sigma_loops_stoch_claude.cu` (l=0, nhits=7, mrhs, at_from_ensdir)
   + `run_sigma_loops_claude.sh` (handoff, user runs).
4. **$F^2$ dump** -- `O` dataset in `glue_f2_v2_shapes_claude.cu` + re-measure handoff.
5. **Coupled GEVP** -- `sigma_f2_gevp_analysis_claude.cu` + run script.
6. **(Later) connected four-point / valence** -- separate.

## Open questions
- **Q-basis:** all 14 $F^2$ shapes in the GEVP, or a reduced set, alongside the 3 pairings + vector$^2$?
- **Q-hits:** nhits=7 the sweet spot, or tune after the gate check tells us the stochastic:gauge balance?
- **Q-scope:** first ensembles (L1/L2, gsq, at0.2)? stride-10 grid confirmed?
- **Q-mrhs width:** block all 7 hits in one solve, or 2 blocks (memory/occupancy)?
