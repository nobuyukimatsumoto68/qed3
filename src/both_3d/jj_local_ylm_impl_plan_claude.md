# Impl plan — $Y_{\ell m}$ tower of the LOCAL vector current correlator (deter ground truth + stochastic conn/disc), per-$m$

Measure the spherical-harmonic descendant tower of the conserved-current two-point function using the
**bare local $\sigma_a$ current** (no conserved-$K$ apply): a deterministic ground truth plus two
stochastic production programs (connected + disconnected), all **per-$m$**, **vector only**, $\ell\le3$.
Designed for the **interacting** ensembles (massless Nf2/4/6 gsq8); the free-field run is the validation
*of that pipeline*.

Reference (CFT tower, validated continuum): `jj_cft_ylm_check_claude.{cc,md}`. Paper: `qed3int_v2-14.pdf`
Eqs. (4.33)-(4.36). Stochastic dilution: Foley et al, hep-lat/0505023. One-end trick: Foster-Michael
(UKQCD), hep-lat/9810021.

## Physics / goal

The tower is the $Y_{\ell m}$ decomposition of the temporal ($\sigma_3$ = $G_t$) and spatial
($(\sigma_1+\sigma_2)/2$ = $G_s$) current correlators between site pairs on the sphere. Diagonal-$\ell$
projection (Funk-Hecke + isotropy, Eq. 4.33):
$$
g_\ell(t)=\tfrac12\int_{-1}^{1} G(t;x)\,P_\ell(x)\,dx,\qquad x=\hat n_1\!\cdot\!\hat n_2 .
$$

CFT prediction (Eq. 4.35, $\Delta=2$, $C_j=-1/(8\pi^2)$ from `jj_cft_ylm_check`):
$$
g^t_1\to\tfrac{1}{24\pi^2}e^{-2t},\qquad g^t_2\to\tfrac{1}{10\pi^2}e^{-3t},\qquad g^t_0\to0 .
$$

Convention-free checks: decay rates $(2,3,4)$ for $\ell=(1,2,3)$; ratio
$[g^t_2 e^{3t}]/[g^t_1 e^{2t}]\to 12/5=2.4$; $g^t_0/g^t_1\to0$; $g^s_\ell/g^t_\ell\to-1$ per $\ell$.

## Estimator (per-$m$, one-end trick)

Per Pauli channel $a\in\{1,2,3\}$ the insertion is the on-site, area-and-harmonic-weighted current
$$
\Sigma^a_{\ell m}(t)=\sum_n A_n\,Y_{\ell m}(\hat n)\,\sigma_a(n,t),
$$
Hermitian; $A_n=$ `dual_areas[n]`, real $Y_{\ell m}$ (`Ylm_real`, `valence_claude.h`); **no $1/\kappa$** (bare
local current). With $\phi=D_m^{-1}\eta$ and per-channel source $\psi^a_{\ell m}=D_m^{-\dagger}\,
\Sigma^a_{\ell m}(t_0)\,\eta$:
$$
g^a_{\ell m}(t)=\psi^{a\,\dagger}_{\ell m}\,\big[\Sigma^a_{\ell m}(t)\,\phi\big]
\quad(\text{stored per-}m),\qquad
g^a_\ell(t)=\frac{1}{2\ell+1}\sum_{m=-\ell}^{\ell} g^a_{\ell m}(t)\ \ (\text{analysis}).
$$
tp tower $=g^3_\ell$; sp tower $=(g^1_\ell+g^2_\ell)/2$. The sink $\Sigma^a_{\ell m}(t)\phi$ is element-wise
(on-site $\sigma_a$ + weight + contraction) — **no Dirac/$K$ apply**, that is the whole cost win over the
conserved current.

## Design decisions (all 2026-06-16, user-confirmed)

- **Single source origin** $t_0$ (default 0). No multiple $t_0$, no master-field superposition; the
  correlator is the ordinary $G(t-t_0)$. Statistics come from configs + hits.
- **Separate conn and disc into two programs** (B1, B2). They need different sources anyway.
- **CONN source = single-timeslice wall** at $t_0$ (`fill_z2_wall_source`); variance-optimal (noise
  concentrated at the origin). $E[\eta\eta^\dagger]=P_{t_0}$, and $\Sigma_{t_0}$ projects onto $t_0$ so the
  connected trace is exact.
- **DISC source = time+spin dilution sweep** (model `disc_claude.cu`): loop `t_s`,`spin` over
  `eta.time_spin_dilution(t_s,t_block,spin)`, `--disc-tblock` default 8. Needed because the one-point
  $J^a_{\ell m}(t)=\mathrm{tr}[\Sigma^a_{\ell m}(t)D_m^{-1}]$ must be defined at every $t$ (the conn wall
  source gives $J$ only at $t_0$). A one-shot full-volume source is unbiased but has ruinous variance —
  dilution is what makes $J(t)$ usable; `t_block` is the variance/cost knob.
- **Store per-$m$**; do the $\frac{1}{2\ell+1}\sum_m$ in analysis. Free in compute (already per-$(a,\ell,m)$).
  Diagnostic: at $L{=}1$ the symmetry is only icosahedral, so $Y_{\ell m}$ modes within an $\ell$ are not
  degenerate — per-$m$ output exposes that anisotropy.
- **$\sigma_{\rm PS}$ condensate folded into the disc program** ($a{=}0$ scalar, uniform-area,
  spacetime-summed loop $=\mathrm{tr}[A D_m^{-1}]$); shares the disc $\phi$, **no extra solve**.
  **$\sigma_{\rm PS}$ only — no $\sigma_{\rm FS}$.**
- **Vector only — no axial.** Massless / $m_F$ target; **$m_P$ (tilde-$D$ backward solve) out of scope**.
- **No block-$t$.** The local sink has no operator to amortize; the $t$-loop is $O(N_t N_{\rm sites})$
  trivial. Cost levers: conn = mrhs over the 16 $(\ell,m)$ columns per Pauli; disc = `t_block`.

## Cost (rough, $\ell\le3$ $\Rightarrow$ 16 $(\ell,m)$, 3 Pauli $\Rightarrow$ 48 channels)

- CONN: per pattern $1$ ($\phi=D_m^{-1}\eta$) $+48$ (source $\psi^a_{\ell m}$) $=49$ solves;
  spin-only (2 patterns) $\approx 98$ solves/hit. mrhs batch (16 cols/Pauli) for ~3-4x.
- DISC: $\tfrac{N_t}{t_block}\times n_\text{spin}\approx 32$ solves/hit at $t_block{=}8$; $\phi$ shared
  across all $(a,\ell,m)$ and the condensate (only the readout weight differs). Disc is intrinsically
  noisy (config-average-dominated); $=0$ in free.

## Chunk A — deterministic ground truth (`jj_local_deter_claude.cu`)

The deterministic per-$m$ ylm is already correct (`build_Sigma_ylm`, $\sigma_3$/s3 only, $\ell\le2$). Two
generalizations:
- Generalize `build_Sigma_ylm` to a **Pauli index `a`** (use `DW.sigma[a]` 2x2 block scaled by
  $A_n Y_{\ell m}$, exactly like `local_W_sigma`); the current $\sigma_3$ hardcode becomes `a==3`. Comment
  the old signature/body, add the generalized one beneath (per convention).
- Tower loop `for a in {1,2,3}` / `for l in 0..L_MAX_YLM` (bump `L_MAX_YLM = 2` -> `3`, comment old) /
  `for m in -l..l`. `Ylm_real` (`std::assoc_legendre`) already valid $\forall\ell$. Store **per-$m$**
  $g^a_{\ell m}=\mathrm{tr}[\Sigma^a_{\ell m}(t_0)P\Sigma^a_{\ell m}(t)P]$ (no $1/(2\ell+1)$, no $m$-sum).
- Output `h0/t0_b/ylm/s{a}/l{l}/m{m}/{Vpp,Vmm}` (supersedes the old s3-only, $m$-summed `ylm/l{l}`).
- Rerun free L=1 ground truth via a `tmp_*_claude.sh` (user runs; existing `corr_deter_local_L1` is
  `complete`-gated -> user `rm`s it first). L=2/L=4 optional later.

## Chunk B1 — stochastic CONN (`jj_local_ylm_conn_stoch_claude.cu`, new)

Solve/output scaffolding from `jj_corr_dilute_claude.cu` (VECTOR path); wall-source + $Y_{\ell m}$/$\sigma$
mechanics from `meson_pq_wall_v2_claude.cu` (`fill_z2_wall_source`, `mult_Ylm_real`, `mult_sigma`).
- Conn-ylm only: no tp/sp/local-s{1,2,3}/axial/disc. Optional spin dilution (`--spin-dilution`); per-hit
  atomic + `complete`-gated output; RNG seed-from-string. Single `--t0` (default 0).
- `ylw[l][m][n] = dual_areas[n]*Ylm_real(l,m,base.sites[n])`, $\ell=0..3$.
- Source $\eta$ = wall at $t_0$; $\phi=D_m^{-1}\eta$. For each `a in {1,2,3}`, each $(\ell,m)$:
  source $\Sigma^a_{\ell m}(t_0)\eta$ (`mult_sigma(a)` + `mult_Ylm_real(l,m)` scaled by `ylw`, restricted to
  $t_0$), solve $\psi^a_{\ell m}=D_m^{-\dagger}(\cdot)$ (`op_Dm`+`op_Dmsq`). Sink loop over $t$:
  $\Phi=\Sigma^a_{\ell m}(t)\phi$; store **per-$m$** `Gv[a,l,m][t-t0] = psi.dag(Phi)`.
- Output `h0/ylm/s{a}/l{0..3}/m{-l..l}/Vpp`.
- (Optional) mrhs-batch the 48 source solves per Pauli via BlockedMat.

## Chunk B2 — stochastic DISC + $\sigma_{\rm PS}$ (`jj_local_ylm_disc_stoch_claude.cu`, new)

The one-point / loop program; ylm-generalization of `saved_scripts_claude/disc_claude.cu` (its loop
structure, `compute_disc_corr`, `--t_block`). Vector disc currents AND $\sigma_{\rm PS}$ share the same
$\phi=D_m^{-1}\eta$ (no extra solve for the condensate).
- `--disc-tblock` (default 8). Loop `t_s=0..Nt/t_block-1`, `spin`:
  `eta.time_spin_dilution(t_s,t_block,spin)`; $\phi=D_m^{-1}\eta$ (per pattern, shared by all readouts).
- DISC currents: per $(a,\ell,m)$, $a\in\{1,2,3\}$, accumulate per-$m$
  $J^a_{\ell m}[t]\mathrel{+}=\sum_n A_n Y_{\ell m}(\hat n)\,\eta^*(t,n,\text{spin})(\sigma_a\phi)
  (t,n,\text{spin})$ at the diluted $t$ — a `ylw`-weighted variant of `accumulate_loop_gamma`.
  Output `h0/disc/ylm/s{a}/l{0..3}/m{-l..l}/J` (RAW per-$m$ $J(t)$, summed over patterns; analysis
  $C_{\rm disc},\ell=\frac{1}{2\ell+1}\sum_m$ `two_point`($J^a_{\ell m}$)).
- CONDENSATE $\sigma_{\rm PS}$: scalar ($a{=}0$, identity), uniform-area, spacetime-summed loop
  $=\mathrm{tr}[A D_m^{-1}]$. Accumulate `etadag_xi += eta_A.dag(phi)` with `eta_A = A*eta` (mirrors the
  dilute program's `acc_etadag_phi`), summed over patterns. $\sigma_{\rm PS}=$ `etadag_xi` $+$ its
  conjugate (backward leg via cyclicity; valid for massless/$m_F$). Output `h0/condensate/etadag_xi`.
  Contact subtraction (`condensate_contact_massive_claude.md` Sec 10) applied in analysis.

## Chunk C — run scripts + validation notebook

Files: `tmp_ylm_conn_stoch_claude.sh`, `tmp_ylm_disc_stoch_claude.sh`, `jj_local_ylm_validate_claude.ipynb`
(all new).
- Build DISTINCT binaries (`jj_local_ylm_conn_stoch.o`, `jj_local_ylm_disc_stoch.o`) to avoid ETXTBSY.
- Free spin-only nhits=140, single $t_0=0$; conn -> `corr_ylm_conn_*`, disc -> `corr_ylm_disc_*`.
- CONN check: load `h0/ylm/s{a}/l{0..3}/m{m}/Vpp`, jackknife over hits; form per-$m$ then
  $g_\ell=\frac{1}{2\ell+1}\sum_m$; tp $=s_3$, sp $=(s_1+s_2)/2$; overlay deterministic
  `corr_deter_local_L1`. Checks: rates $\to(2,3,4)$; $g^t_2 e^{3t}/g^t_1 e^{2t}\to2.4$; $g_0\to0$;
  sp/tp $\to-1$; **stoch == determ** within jk error (per-$m$ AND $\ell$-summed). Honest normalization
  (one constant per estimator, raw ratios, no np.abs / per-curve sign flips).
- DISC check: $C_{\rm disc}=$ `two_point`($J$) $\to0$ within jk error in free; physical $=-C_{\rm conn}+
  C_{\rm disc}$ (Eq. 3.39 sign). Mirrors the disc section in the `jj_validate_*` notebooks.
- CONDENSATE check: $\sigma_{\rm PS}\to-2$ in free after contact subtraction (cf. dense
  `condensate_deter_free`).

## Reference files

- `jj_local_deter_claude.cu` — existing per-$m$ deterministic ylm (the always-correct reference); Chunk A
  edits this.
- `meson_pq_wall_v2_claude.cu` — wall source + `mult_Ylm_real`/`mult_sigma` + one-solve-per-channel
  (conn template).
- `saved_scripts_claude/disc_claude.cu` — volume time+spin dilution loop + `accumulate_loop_gamma` +
  `compute_disc_corr` (disc template).
- `jj_corr_dilute_claude.cu` — solve operators (`op_Dm`/`op_Dmsq`), output/atomic/complete-gate, condensate
  `acc_etadag_phi` (scaffolding).
- `jj_cft_ylm_check_claude.{cc,md}` — CFT tower targets.

## Open questions

- **$\ell=3$ amplitude.** Validate rate$\to4$ and stoch==determ; a new CFT amplitude constant for
  $g_3/g_2$ is optional (extend `jj_cft_ylm_check` to $\ell=3$ only if wanted).
- **`t_block` for disc** in the interacting run — default 8 (already tuned on Nf2 gsq8 in
  `disc_Nf2_gsq8...tb8/`); revisit if the disc signal is too noisy.

## Addendum (2026-06-16) — AXIAL Y_lm tower (no GW dressing)

User: add an axial tower alongside the vector one, modeled on the LOCAL axial in `jj_corr_dilute_claude.cu`
(`:778-807`) — **NOT** the exact-K connected axial (which has the $(1-D_{ov})$ GW factor). The local axial
is the local VECTOR with the $t_0$-leg propagator swapped from backward to forward, bare $\sigma_a$ vertex:

- dilute local vector source (`:743`): `op_Dm` RHS -> $\psi = D_m^{-\dagger}(\sigma_c\eta)$ (backward leg).
- dilute local axial source (`:794`): `op_DmH` RHS -> $\psi_A = D_m^{-1}(\sigma_c\eta)$ (forward leg).
- SAME sink $\chi=\sigma_c\phi'$, $\phi'=D_m^{-1}\eta$.  NO $(1-D_{ov})$.

Ylm axial estimator (per $a,\ell,m$; same $\Sigma^a_{\ell m}$, $\phi=D_m^{-1}\eta$ as vector):
$$
\psi^a_{\ell m,A}=D_m^{-1}\,\Sigma^a_{\ell m}(t_0)\,\eta,\qquad
g^a_{\ell m,A}(t)=\psi^{a\,\dagger}_{\ell m,A}\,[\Sigma^a_{\ell m}(t)\,\phi]
=\mathrm{tr}[\Sigma_0\,P^\dagger\Sigma_t\,P]\quad(P=D_m^{-1}).
$$
One dagger vs the vector's $\mathrm{tr}[\Sigma_0 P\Sigma_t P]$.  Determ: a `tr_WPdagWP` variant indexing
`conj(P[et.i,e0.j])` on the first leg (no extra memory).  m_P (parity) still OUT OF SCOPE.

**Implementation:**
- `jj_local_deter_claude.cu` — add `tr_WPdagWP`; axial tower loop -> `h0/t0_b/ylm_axial/s{a}/l{l}/m{m}/{Vpp,Vmm}`.
- `jj_local_ylm_conn_stoch_claude.cu` — per $(a,\ell,m)$ also solve $\psi_A$ (forward: `op_DmH`+`op_Dmsq`,
  same as $\phi$) and accumulate -> `h0/ylm_axial/s{a}/l{l}/m{m}/Vpp`.  +48 source solves/pattern (~194/hit).
- **DISC: none for axial** — the one-point $\mathrm{tr}[\sigma_c S]$ is vertex-only (= vector); the local
  axial is connected-only in the reference.  B2 unchanged.
- Notebook: axial tower cells (signed-log, stoch vs determ), parallel to the vector tp/s1/s2 cells.
