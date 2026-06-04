# Plan v2: Conserved-Current Correlator Measurements

Clean re-summary of the design worked out in
`conserved_current_correlators_impl_plan_v1_claude.md` (v1, scratch).  All decisions below
are settled; v1 keeps the derivation history.

Reference: A. Katz, N. Matsumoto, "Notes on QED3 (interacting theory)",
`/mnt/barracuda22/qed3/qed3int_v2-10.pdf` (renewed; supersedes `_v2-9`, `_v2`) --
Secs. 3.3.1--3.3.4 (lattice correlators, Eqs. 3.39--3.69), 4.2--4.3 (CFT estimators,
Eqs. 4.27--4.36; Section-4 numbering unchanged), App. C (real spherical harmonics).

---

## 1. Goal and scope

Measure current-current two-point functions on the $S^2\times$time lattice using the
exactly conserved overlap currents, for $N_f = 2, 4, 6$.  Currents are two-component in the
$\xi$/$\eta$ basis; the (+) vector current is $J_{V,+}^{nn'}(t) = \eta_f^\dagger K^{nn'}(t)\xi_f$.

Current types (PDF naming, Sec. 3.3.1, Eqs. 3.39--3.44):

- **Vector** $J_{V,\pm}$ -- connected $+$ disconnected.  The $(+)$ correlator gives
  $C_{V_{++},c}$, $C_{V_{++},d}$; the $(-)$ gives $C_{V_{--},c}$, $C_{V_{--},d}$ (massless and
  $m_F$: $C_{V_{--}}=C_{V_{++}}^*$, Eq. 3.54; parity $m_P$: independent, Eqs. 3.66--3.67).
- **Flavor** $J_{r,\pm}$ -- connected only, equal to the vector-connected piece
  $C_{V_{\pm\pm},c}$ (Eqs. 3.43--3.44); **does not exist for $N_f = 2$** (relevant for
  $N_f = 4, 6$).
- **Axial** $J_{A,\pm}$ -- connected only, the two *mixed* orderings
  $C_{A_{+-}}(0\to t)\propto\langle J_{A,+}(0)J_{A,-}(t)\rangle_F$ and
  $C_{A_{-+}}(0\to t)\propto\langle J_{A,-}(0)J_{A,+}(t)\rangle_F$ (Eqs. 3.41--3.42), related
  by $C_{A_{+-}}(0\to t)=C_{A_{-+}}(t\to 0)$ (Eq. 3.57).

Three projected observables, computed for each current type:

- $G_s(t)$ -- spatial, equal-link average (Eq. 4.29).
- $G_t(t)$ -- temporal, site-diagonal (Eq. 4.32).
- $G^{t}_{\ell_1 m_1;\ell_2 m_2}(t)$ -- spherical-harmonic projection (Eq. 4.36), which
  carries the conformal tower (Eqs. 4.34--4.35): the $\ell=0$ channel vanishes (charge
  conservation), $\ell=1 \sim e^{-2t}$, $\ell=2 \sim e^{-3t}$.

Mass cases (all supported by every program, PDF Secs. 3.3.2--3.3.4):

- **massless** ($m_F=m_P=0$; Eqs. 3.50--3.57);
- **flavor-breaking** $m_F\in\mathbb{R}$ (real, parity-symmetric; Eqs. 3.58--3.60);
- **parity-breaking** $m_P\in i\mathbb{R}$ (purely imaginary, flavor-symmetric;
  Eqs. 3.61--3.69).

The propagator operators for each case ($D_m=D_\text{ov}+m$, and the dagger-leg
$\tilde D_{m_P}=D_\text{ov}+m_P/(1-m_P)$ with general complex mass) are in Sec. 3.7; both are
in the $D_\text{ov}+m$ form, so the existing solver routines apply unchanged.

### 1.1 Per-ensemble analysis matrix

The ensembles live on a remote cluster, several at real ($m_F$) and imaginary ($m_P$) sea
masses.  "Sea" = the mass the gauge configs were generated with (selects which configs to
read); "valence" = the mass passed to the measurement executable (drives all propagators).
They are decoupled and differ only in the $m_F$-axial cell:

| Ensemble (sea) | Vector (conn + disc) | Axial (conn only; no disc) |
|---|---|---|
| massless | valence $0$ | valence $0$ |
| $m_F$ (real) | valence $m_F$ -- exact massive current (Eqs. 3.58--3.60) | valence $0$ -- massless formulas (Eqs. 3.55--3.57) on the $m_F$ configs |
| $m_P$ (imag) | valence $m_P$ -- exact massive current (Eqs. 3.64--3.67) | valence $m_P$ -- exact massive current (Eqs. 3.68--3.69) |

Only the **valence propagator mass** varies across cells; the kernel $K$ is always the same
massless-form Noether kernel (there is no "kernel mass").  Any finite correction to $K$ is
exactly compensated by a propagator factor (the $(1+m_P)^{-1}$ of Eq. 3.62), so the
implemented traces use the standard $K$ with **no prefactor and no case branching** in the
executables.  Disconnected is **vector-only** (axial has no disconnected diagram).  At finite
$m_P$ all continuous symmetries stay intact, so the axial current is still exactly conserved
(unlike $m_F$, where no exactly conserved axial current exists -- hence the massless form).

A driver **shell script** (Chunk 10) holds the case logic: it iterates ensembles
(gsq, $N_f$, sea mass) and, per row above, points each executable at the sea config directory
while passing the appropriate valence `--mass-re`/`--mass-im` (and `--current`).

---

## 2. Notation

Fermionic contractions (Eqs. 3.50--3.51), with the fermionic average $\langle\cdot\rangle_F$:

$$
\langle \xi_{f,x_1}\,\eta_{g,x_2}^\dagger\rangle_F = \delta_{fg}(D_\text{ov}^{-1})_{x_1 x_2},
\qquad
\langle \eta_{f,x_1}\,\xi_{g,x_2}^\dagger\rangle_F = \delta_{fg}(D_\text{ov}^{-\dagger})_{x_1 x_2}.
$$

There is **no $\gamma_5$-hermiticity** here ($D_\text{ov}^{-\dagger}\neq\gamma_5 D_\text{ov}^{-1}\gamma_5$).
$\langle\cdot\rangle_U$ is the gauge average; $^*$ is complex conjugation.  The two contractions
above are the massless ($m_F=m_P=0$) case; in the massive cases $D_\text{ov}$ and
$D_\text{ov}^\dagger$ are replaced by the case-dependent operators of Sec. 3.7, and the
substitution propagates through every estimator of Secs. 3.1--3.4.

Kernels: $K^{nn'}(t)$ (spatial link) and $K^{t,t+1}(n)$ (temporal link), via
`ConservedCurrent::apply_k`; the adjoint $K^\dagger$ via `apply_k_dag`
(`conserved_current_claude.h:157`, spatial wrapper `apply_k_dag_wz` `:230`; implemented,
under verification).  Both accept `std::pair<int,BaseLink>` (spatial) and
`std::pair<int,Idx>` (temporal).

---

## 3. Estimators

### 3.1 Disconnected (vector only) -- factorizes

$$
C_{V_{++},d}^{nn'}(0\to t) = T^{nn'}(0)\,T^{nn'}(t),
\qquad
T^{nn'}(t) \equiv \mathrm{tr}[D_\text{ov}^{-1}K^{nn'}(t)] \approx \eta^\dagger K^{nn'}(t)\,\phi,
\quad \phi = D_\text{ov}^{-1}\eta,
$$

(Eq. 3.53) estimated with $Z_2\times Z_2$ noise and time-spin dilution (as `disc_claude.cu`).
The $(-)$ channel is $C_{V_{--},d}=C_{V_{++},d}^*$ (massless/$m_F$, Eq. 3.54), so the same
$T$ suffices.  **Parity caveat:** in the $m_P$ case the $(-)$ disconnected uses the
*dagger-leg* trace, $C_{V_{--},d}=\tilde T(0)\tilde T(t)$ with
$\tilde T\equiv\mathrm{tr}[\tilde D_{m_P}^{\dagger-1}K^\dagger]$ (Eq. 3.67, no prefactor in
v2-10), so the dump must also emit $\tilde T$ then (Sec. 3.7).

The measurement program **dumps the raw per-config single-time traces** $T^{nn'}(t)$ (spatial)
and $T^{t,t+1}(n)$ (temporal); the product, gauge average, vacuum subtraction, weighting,
and projection are done later (Secs. 3.4--3.5, 5).

### 3.2 Connected (vector / flavor) -- single source, source/sink split

$$
\phi = K^{(\text{src})}\,D_\text{ov}^{-1}\eta,
\qquad
\psi = D_\text{ov}^{-\dagger}\,K^{(\text{snk})\dagger}\,\eta,
\qquad
C_{V_{++},c} = \langle \psi^\dagger\phi\rangle .
$$

Since $\psi^\dagger = \eta^\dagger K^{(\text{snk})}D_\text{ov}^{-1}$,
$E[\psi^\dagger\phi] = \mathrm{tr}[D_\text{ov}^{-1}K^{(\text{snk})}D_\text{ov}^{-1}K^{(\text{src})}] = C_{V_{++},c}$
(Eq. 3.52): a single bilinear, exact, both legs forward $D_\text{ov}^{-1}$, no disconnected
contamination.  The $(-)$ channel is $C_{V_{--},c}=C_{V_{++},c}^*$ (massless/$m_F$, Eq. 3.54);
parity is Eq. 3.66 (no prefactor in v2-10) with the dagger-leg operators of Sec. 3.7.

Efficiency: the **source side is free** -- $\phi_0 = D_\text{ov}^{-1}\eta$ is solved once
per hit and the source kernel is applied after it, covering all source insertion points and
times.  The **sink side** is one $D_\text{ov}^\dagger$ solve per sink time on a *composite*
source (Sec. 3.4).

### 3.3 Axial (mixed orderings) -- same split with GW factors

$C_{A_{+-}}(0\to t)=\mathrm{tr}[D_\text{ov}^{-1}K(0)(1-D_\text{ov}^\dagger)D_\text{ov}^{-\dagger}K^\dagger(t)(1-D_\text{ov})]$
(Eq. 3.55):

$$
\phi_A = D_\text{ov}^{-\dagger}\,K^\dagger(t)\,(1-D_\text{ov})\,\eta \ (\text{sink side}),
\qquad
\psi_A = (1-D_\text{ov})\,K^\dagger(0)\,\rho_0 \ (\text{source side, free}),
\qquad \rho_0 = D_\text{ov}^{-\dagger}\eta,
$$

with $C_{A_{+-}}(0\to t) = \langle\psi_A^\dagger\phi_A\rangle$.  The other ordering
$C_{A_{-+}}(0\to t)$ (Eq. 3.56) is the $K\leftrightarrow K^\dagger$,
$D_\text{ov}\leftrightarrow D_\text{ov}^\dagger$ mirror (source side uses
$\phi_0=D_\text{ov}^{-1}\eta$); equivalently $C_{A_{+-}}(0\to t)=C_{A_{-+}}(t\to 0)$
(Eq. 3.57).  The GW factors $(1-D_\text{ov})$, $(1-D_\text{ov}^\dagger)$ stay **massless**
even for $m_F$/$m_P$ (parity has no explicit prefactor in v2-10, Eqs. 3.68--3.69).
Implemented **as written**.

### 3.4 Composite sink source (one solve per sink time)

The sink-side solve always acts on a composite source with the insertion-point sum folded
in, so it is **one $D_\text{ov}^\dagger$ solve per sink time**, never per insertion point:

- **`_ylmproj`** -- fold the deterministic $Y_{\ell m}$ weights:
  $\Psi^{(\ell_2 m_2)}(s) = D_\text{ov}^{-\dagger}\sum_{n}\frac{A_{n}Y_{\ell_2 m_2}(\hat n)}{\kappa^{(0)}_{t,t+1}(n)}K^{t,t+1\,\dagger}(n,s)\eta$
  (one solve per channel per sink time), paired with
  $\Phi^{(\ell_1 m_1)}(s_0)=\sum_n\frac{A_n Y_{\ell_1 m_1}(\hat n)}{\kappa^{(0)}_{t,t+1}(n)}K^{t,t+1}(n,s_0)\phi_0$ (free).
- **`_spproj` / `_tpproj` (diagonal)** -- realize the diagonal sum with stochastic
  insertion-point noise $\zeta_a$ ($a$ = link / site, $E[\zeta_a\zeta_b^*]=\delta_{ab}$):

$$
\Phi^\zeta(s_0)=\sum_a \zeta_a\sqrt{w_a}\,K^a(s_0)\,\phi_0 \ (\text{free}),
\qquad
\Psi^\zeta(s)=D_\text{ov}^{-\dagger}\sum_a \zeta_a\sqrt{w_a}\,K^{a\dagger}(s)\,\eta \ (\text{one solve}),
$$

  giving $E_{\zeta,\eta}[\Psi^\zeta(s)^\dagger\Phi^\zeta(s_0)] = \sum_a w_a\,\mathrm{tr}[K^a(s)D_\text{ov}^{-1}K^a(s_0)D_\text{ov}^{-1}]$.
  Off-diagonal terms average to zero; their variance is controlled by the number of
  $\zeta$-hits.

Because $\Phi$ supplies all source times for free, solving $\Psi$ at $n_s$ sink times yields
all $\Delta t$ with $n_s$ source-time samples.

### 3.5 Section 4 observables (geometric weights applied in-program)

$$
G_s(t) = \frac{1}{4\pi}\sum_{\langle n,n'\rangle}\frac{A_{nn'}}{\kappa^{(0)\,2}_{nn'}}\langle J^{nn'}(0)J^{nn'}(t)\rangle,
\qquad
G_t(t) = \frac{1}{4\pi}\sum_{n}\frac{A_{n}}{\kappa^{(0)\,2}_{t,t+1}(n)}\langle J^{0,1}(n)J^{t,t+1}(n)\rangle,
$$

$$
G^{t}_{\ell_1 m_1;\ell_2 m_2}(t) = \frac{1}{4\pi}\sum_{n_1,n_2}
  \frac{A_{n_1}Y_{\ell_1 m_1}(\hat n_1)}{\kappa^{(0)}_{t,t+1}(n_1)}
  \frac{A_{n_2}Y_{\ell_2 m_2}(\hat n_2)}{\kappa^{(0)}_{t,t+1}(n_2)}
  \langle J^{0,1}(n_1)J^{t,t+1}(n_2)\rangle .
$$

For the vector $(+)$ current $\langle JJ\rangle = -C_{V_{++},c}+C_{V_{++},d}$ (and the
$(-)$ analogue).  The geometric weights (and the real $Y_{\ell m}$ from `GetYlm`) are applied
**inside the compiled programs**; the notebooks only gauge-average, vacuum-subtract, combine,
and fit.

**Normalization convention.**  The $C$'s above are the bare traces; the PDF relations carry
overall factors, e.g. $\tfrac{2}{N_f}\langle J_{V,+}(0)J_{V,+}(t)\rangle_F = -C_{V_{++},c}+C_{V_{++},d}$
(Eq. 3.39) and the disconnected loop is $L=-N_f\,T$.  These overall $2/N_f$, $-N_f$ factors
are bookkeeping applied at the analysis stage; the programs output the bare $C_{V_{\pm\pm},c}$,
$C_{A_{\pm\mp}}$, $T$.

### 3.6 Geometric quantities (already in code)

| Symbol | Meaning | Accessor |
|---|---|---|
| $\kappa^{(0)}_{nn'}$ | free-field spatial-link coupling | `DiracOp::kappa[il]` (`dirac_simp.h:355-360`) |
| $A_{nn'}$ | dual area of edge $(n,n')$ | `lattice.link_volume[il]` (`s2n_simp.h:29`) |
| $A_n$ | dual (Voronoi) area of site $n$ | `lattice.dual_areas[s]` (`s2n_simp.h:38`) |
| $\kappa^{(0)}_{t,t+1}(n)$ | temporal-link coupling (uniform in $t$, per site) | per-site array (candidate `dual_areas[ix]`); pin when coding |
| $Y_{\ell m}(\hat n)$ | real spherical harmonic | `lattice.GetYlm(s,l,m)` after `UpdateYlm(l_max)` (`s2.h:55-56`) |

### 3.7 Mass cases: massless / flavor-breaking $m_F$ / parity-breaking $m_P$

The programs support three mass cases (PDF Secs. 3.3.2--3.3.4).  Define (Eqs. 3.60, 3.63)

$$
D_m \equiv D_\text{ov}+m,
\qquad
\tilde D_{m_P} \equiv D_\text{ov}+\frac{m_P}{1-m_P}.
$$

Both have **coefficient 1 on $D_\text{ov}$**, i.e. the $D_\text{ov}+m$ form.  In every
estimator of Secs. 3.1--3.4 the symbols $D_\text{ov}$, $D_\text{ov}^\dagger$ are replaced
case-by-case as below; the kernels $K$, $K^\dagger$ and the geometric weights are unchanged.

> **Definition swap vs old memory.**  As in the PDF, $m_F$ is the **real**, parity-symmetric,
> *flavor-breaking* mass and $m_P$ is the **purely imaginary**, flavor-symmetric,
> *parity-breaking* mass -- **swapped** relative to the older memory notes.

**Massless** ($m_F=m_P=0$; Eqs. 3.50--3.57).
$\langle\xi\eta^\dagger\rangle_F=D_\text{ov}^{-1}$,
$\langle\eta\xi^\dagger\rangle_F=D_\text{ov}^{-\dagger}$.

**Flavor-breaking** ($m_F\in\mathbb{R}$, $m_P=0$; Eqs. 3.58--3.60).

$$
\langle\xi\eta^\dagger\rangle_F=D_{m_F}^{-1},
\qquad
\langle\eta\xi^\dagger\rangle_F=(D_{m_F}^{-1})^\dagger .
$$

Substitute $D_\text{ov}\to D_{m_F}$ throughout Eqs. 3.52--3.57.  Since $m_F$ is real,
$(D_{m_F})^\dagger=D_\text{ov}^\dagger+m_F$ and the dagger leg is automatically consistent.

**Parity-breaking** ($m_P\in i\mathbb{R}$, $m_F=0$; Eqs. 3.61--3.69).

$$
\langle\xi\eta^\dagger\rangle_F=D_{m_P}^{-1},
\qquad
\langle\eta\xi^\dagger\rangle_F=(1+m_P)^{-1}(\tilde D_{m_P}^{-1})^\dagger ,
\qquad m'\equiv\frac{m_P}{1-m_P}.
$$

The forward leg uses $D_{m_P}=D_\text{ov}+m_P$ (`mass`$=m_P$, purely imaginary); the dagger
leg uses $\tilde D_{m_P}=D_\text{ov}+m'$ (`mass`$=m'=m_P/(1-m_P)$, **general complex**:
neither real nor purely imaginary, since $m_P=i\mu\Rightarrow m'=(-\mu^2+i\mu)/(1+\mu^2)$).
Because the coefficient of $D_\text{ov}$ is 1, the GW relation holds and the existing
`DHD`/`DDH` routines are usable as-is (see below) -- this is the v2-10 simplification that
replaces the earlier $(1-m_P)D_\text{ov}+m_P$ form.

In the correlators (Eqs. 3.64--3.69) the $(1+m_P)^{-1}$ propagator prefactor **cancels the
finite correction of the adjoint kernel**, so the implemented traces carry **no explicit
prefactor**: $C_{V_{++},c/d}$ use $D_{m_P}^{-1}$; $C_{V_{--},c/d}$ use
$\tilde D_{m_P}^{\dagger-1}K^\dagger$; the axial $C_{A_{+-}},C_{A_{-+}}$ mix $D_{m_P}^{-1}$ and
$\tilde D_{m_P}^{\dagger-1}$ with the GW factors kept **massless** $(1-D_\text{ov})$,
$(1-D_\text{ov}^\dagger)$.

**Solver and infrastructure (`overlap_wmass_claude.h`).**  $D_m$ and $\tilde D_{m_P}$ are
both of the form $D_\text{ov}+m$, so the existing complex-mass overlap (used by
`hmc_w_mass_claude.cu` as `Overlap D(DW, mass, 21)`, `mass=Complex(mass_re,mass_im)`) covers
them with no operator change: `mult`$=(D_\text{ov}+m)v$ (`:334`),
`adj`$=(D_\text{ov}+m)^\dagger v=D_\text{ov}^\dagger v+m^* v$ (`:366`),
`DHD`$=(D+m)^\dagger(D+m)$ (`:371`), `DDH` aliased to `DHD` (`:392`).

- `DHD` already has a **general complex $m$ in scope** (verified): it uses the identity
  $(D{+}m)^\dagger(D{+}m)=(1{+}m^*)(D{+}m)+(1{+}m)(D{+}m)^\dagger-(2\,\mathrm{Re}\,m+|m|^2)$ with
  `std::conj(mass)` (`:381`), `Complex(1.0)+mass` (`:383`), `2.0*mass.real()+std::norm(mass)`
  (`:385`) -- all correct for any complex $m$.  The identity holds for any complex $m$; it
  needs only $D_\text{ov}^\dagger D_\text{ov}=D_\text{ov}+D_\text{ov}^\dagger$ (GW, coefficient
  of $V=D_\text{ov}-1$ equal to 1), which $\tilde D_{m_P}$ now satisfies.
- `DDH=DHD` (`:392`) is valid here because $D_\text{ov}+m$ is normal (a scalar shift of the
  normal $D_\text{ov}$), so $DD^\dagger=D^\dagger D$.  Kept as-is per the v2-10 design.
- A program needing both legs instantiates two `Overlap` objects (or switches `mass`): one at
  `mass`$=m_P$ (forward), one at `mass`$=m'$ (dagger); both share `DW` and `update(U)`.

---

## 4. Programs and data flow

| File | Role |
|---|---|
| `jj_disc_claude.cu` | **GPU dump**: per-config raw traces $T^{nn'}(t)$ (spatial) and $T^{t,t+1}(n)$ (temporal).  One program. |
| `jj_conn_spproj_claude.cu` | **GPU** connected, Eq. 4.29 (spatial).  `--current vector\|axial`. |
| `jj_conn_tpproj_claude.cu` | **GPU** connected, Eq. 4.32 (temporal diagonal).  `--current`. |
| `jj_conn_ylmproj_claude.cu` | **GPU** connected, Eq. 4.36 ($Y_{\ell m}$).  `--current`. |
| `jj_disc_spproj_claude.cu` | **CPU/GPU** disc post-proc, Eq. 4.29 (reads the dump, applies weights). |
| `jj_disc_tpproj_claude.cu` | disc post-proc, Eq. 4.32. |
| `jj_disc_ylmproj_claude.cu` | disc post-proc, Eq. 4.36. |
| `analysis/*.ipynb` | gauge average + jackknife, vacuum subtraction, $-C_{V,c}+C_{V,d}$ combination, fits/plots. |

Data flow: GPU dump ($T$) $\to$ disc post-proc ($G^{disc}$ per config, weighted) ; GPU
connected ($G^{conn}$ per config, weighted, per `--current`) ; both $\to$ notebooks.

Conventions:

- All geometric weighting and $Y_{\ell m}$ projection happen in the **compiled programs**;
  notebooks do only averaging / vacuum subtraction / fitting.
- `--current vector|axial` selects the current for the connected programs (one per run; no
  inversion reuse across currents -- simplicity).
- **Valence mass** (all programs): complex `--mass-re`/`--mass-im` = the propagator mass
  ($0\Rightarrow$ massless, real$\Rightarrow m_F$, imaginary$\Rightarrow m_P$).  Forward solve
  uses $D_m=D_\text{ov}+m$; the parity dagger leg uses $\tilde D_{m_P}=D_\text{ov}+m_P/(1-m_P)$
  (coefficient 1) via the standard DHD/DDH solve (Sec. 3.7) -- existing routines, general
  complex mass.  There is **no kernel flag** and no case branching: $K$ is always the
  massless-form kernel (Sec. 1.1).
- **Ensemble (sea) selection**, decoupled from the valence mass: the config directory (or the
  sea-mass identifier that builds it) is its own CLI argument, so a massless-valence run can
  read an $m_F$ ensemble (the $m_F$-axial cell of Sec. 1.1).  Output directory/filename
  encodes **both** the sea mass (ensemble) and the valence mass + current + projection
  (`mRe`/`mIm` as `hmc_w_mass_claude.cu:186`), to avoid collisions when valence $\neq$ sea.
- Production noise: `nhits=1`, `t_block=8` (time-spin dilution), as `disc_claude.cu:164-165`.
  Connected programs use all-to-all $\eta$ (no `--t_block`); diagonal projections add
  $\zeta$ insertion-point noise.
- Output: one HDF5 file per config $k$, `HighFive::File::Truncate`, with resume sentinel
  (as `disc_claude.cu:272-286`).

Compile-time (inherit `disc_claude.cu`): `N_REFINE=1`, `Nt=128`, `NPARALLEL_DUPDATE=1`,
`TOL_INNER=1e-15`, `TOL_OUTER=1e-14`.

Reused infrastructure: complex-mass overlap `Overlap<WilsonDirac>` from
`overlap_wmass_claude.h` (built with `Complex mass`; reduces to massless at `mass=0`) -- use
this in place of the massless `overlap.h` so $D_m$/$\tilde D_{m_P}$ come for free;
`apply_k` / `apply_k_dag` (`conserved_current_claude.h`); CG forward chain `op_DH` /
`op_DHD` (`disc_claude.cu:231-237`), and for the dagger leg a $D_\text{ov}^\dagger$ solve via
the `DDH` chain (`overlap_wmass_claude.h:392`, `DDH=DHD` valid since $D_\text{ov}+m$ is
normal); `from_cpu` / `solve`; `FermionVector` with `fill_z2_source`,
`time_spin_dilution`, `dag` (`valence_claude.h`); `GaugeExt::read` (`disc_claude.cu:282`);
lattice geometry (Sec. 3.6).

---

## 5. Implementation chunks

### Disconnected dump

**Chunk 1 -- boilerplate + construction.**  Copy `disc_claude.cu` boilerplate (includes,
`Comp`, type aliases, `ParseArgs`, GPU/lattice/DW/Gauge setup); use the complex-mass overlap
`overlap_wmass_claude.h` with `Overlap D(DW, mass, 21)`, `mass=Complex(mass_re,mass_im)` and
add `--mass-re`/`--mass-im` to `ParseArgs` (Sec. 3.7); add
`#include "includes/conserved_current_claude.h"` and `ConservedCurrent<Fermion> kop(D)`.
Files: `jj_disc_claude.cu`.

**Chunk 2 -- enumerate + dump layout.**  Enumerate spatial links $(s, lk)$ over `base.links`
and temporal sites $(n, s)$; HDF5 keys `spatial/{ix0}/{ix1}/{real,imag}` (len $N_t$) and
`temporal/{n}/{real,imag}` (len $N_t$); output dir + resume sentinel.
Files: `jj_disc_claude.cu`.

**Chunk 3 -- stochastic loop + dump.**  Per config, per hit, per dilution block: draw
diluted $\eta$, solve $\phi=D_m^{-1}\eta$, accumulate
$T^{lk}(s)\mathrel{+}=\eta^\dagger K^{lk}(s)\phi$ and $T^{n}(s)\mathrel{+}=\eta^\dagger K^{t,t+1}(n,s)\phi$;
divide by `nhits`; write raw $T$ to HDF5.  No correlator formed.  **Parity case only** (Sec.
3.1): also solve the dagger leg $\tilde\phi=(\tilde D_{m_P}^{-1})^\dagger\eta$ and dump
$\tilde T\mathrel{+}=\eta^\dagger K^\dagger\tilde\phi$ for the $(-)$ disconnected (Eq. 3.67).
Files: `jj_disc_claude.cu`.

### Connected (GPU)

**Chunk 4 -- boilerplate + `--current`.**  As Chunk 1 for the three `jj_conn_*proj` programs;
all-to-all $\eta$; add `--current vector|axial` (branch in the loop).
Files: `jj_conn_spproj_claude.cu`, `jj_conn_tpproj_claude.cu`, `jj_conn_ylmproj_claude.cu`.

**Chunk 5 -- per-config loop.**  Per hit: base solve(s) $\phi_0=D_\text{ov}^{-1}\eta$ (and
$\rho_0=D_\text{ov}^{-\dagger}\eta$ for axial).  Build the free source-side vectors (apply
source kernels to $\phi_0$/$\rho_0$); build the composite sink source (Sec. 3.4: $Y_{\ell m}$
fold for ylm, $\zeta$ noise for sp/tp) and do one $D_\text{ov}^\dagger$ solve per sink time;
contract $\psi^\dagger\phi$ over source times to fill $G(\Delta t)$ for the selected current
(vector: $C_{V_{++},c}$, with $C_{V_{--},c}=C_{V_{++},c}^*$ in massless/$m_F$ but a separate
$\tilde D_{m_P}^{\dagger-1}$ trace in parity; axial: $C_{A_{+-}}$ and $C_{A_{-+}}$).  Average
over $\eta$ (and $\zeta$) hits.
Files: the three `jj_conn_*proj_claude.cu`.

**Chunk 6 -- output.**  HDF5 per config; filename encodes `--current` and the mass; datasets:
sp/tp -> $G(t)$ length $N_t$; ylm -> per $(ch_1,ch_2)$ length $N_t$; vector writes
$C_{V_{++},c}$ (and $C_{V_{--},c}$ in parity); axial writes $C_{A_{+-}}$ and $C_{A_{-+}}$.
Resume sentinel.
Files: the three `jj_conn_*proj_claude.cu`.

### Disconnected post-processing

**Chunk 7 -- geometry + read.**  Reconstruct the lattice (Sec. 3.6), build the weight arrays
$w^{sp}_{lk}$, $w^{tp}_n$, $w^{ylm}_{n,(\ell m)}$; enumerate and read `jj_disc_trace.{k}.h5`.
Files: `jj_disc_spproj_claude.cu`, `jj_disc_tpproj_claude.cu`, `jj_disc_ylmproj_claude.cu`.

**Chunk 8 -- projection + output.**  Per config, apply weights in-program and form the
projected products (source-time averaged): sp/tp $\to$ $w$-weighted, summed
$\sum_a w_a\frac1{N_t}\sum_{s_0}T^a(s_0)T^a(s_0+t)$ with per-link/site weighted singles kept
for the vacuum term; ylm $\to$ $\tilde T^{(\ell m)}(s)=\sum_n w^{ylm}_{n,(\ell m)}T^{t,t+1}(n,s)$,
then $D^{ylm}_{(\ell_1 m_1),(\ell_2 m_2)}(t)=\frac1{N_t}\sum_{s_0}\tilde T^{(\ell_1 m_1)}(s_0)(\tilde T^{(\ell_2 m_2)}(s_0+t))^*$.
Write per config in the connected-like layout.
Files: the three `jj_disc_*proj_claude.cu`.

### Validation

**Chunk 9 -- guarded self-checks** (`if constexpr (N_REFINE==1 && Nt<=4)`, cf.
`check_conserved_current_claude.cu:643`): disc -- stochastic $T$ vs exact basis-vector
trace; connected -- $C_{V,c}(\Delta t=0)\geq 0$ and (optional) $\langle\psi^\dagger\phi\rangle$
vs exact $\mathrm{tr}[D^{-1}KD^{-1}K]$.
Files: guarded blocks in the GPU programs (no new file).

### Driver

**Chunk 10 -- ensemble-sweep shell script.**  Driver that runs the analysis executables over
the remote-cluster ensembles, holding the Sec. 1.1 case logic.  Inputs: a list of ensembles
(gsq, $N_f$, sea mass) and their config directories.  Per ensemble it dispatches, following
the matrix: **vector** -> disc dump (`jj_disc`) + the three `jj_conn_*proj` with
`--current vector`; **axial** -> the three `jj_conn_*proj` with `--current axial` (no disc).
It sets the **valence** `--mass-re`/`--mass-im` per cell -- equal to the sea mass for
vector(all) and for $m_P$ axial, but **0** for $m_F$ axial (massless formulas on the $m_F$
configs) -- while always pointing at the sea config directory.  No `if`-branching inside the
executables; all case logic lives here.  Honor the per-program resume sentinels so reruns
skip completed configs.  Sized for the user's cluster workload (not the 4-core cap).
Files: `run_jj_analysis_claude.sh`.

---

## 6. Pending / future

**Open questions (mass cases).**

- **Q1 -- CLI surface.**  RESOLVED: complex `--mass-re`/`--mass-im` matching
  `hmc_w_mass_claude.cu` ($0\Rightarrow$ massless, real$\Rightarrow m_F$,
  imaginary$\Rightarrow m_P$).
- **Q2 -- dagger-side normal operator (DDH).**  RESOLVED by v2-10's coefficient-1
  $\tilde D_{m_P}=D_\text{ov}+m_P/(1-m_P)$ (Eq. 3.63): the dagger operator is in the
  $D_\text{ov}+m$ form, so the existing `DHD`/`DDH` routines are usable unchanged.  `DHD`
  already handles a general complex $m$ (verified, `overlap_wmass_claude.h:381,383,385`), and
  `DDH=DHD` (`:392`) is valid since $D_\text{ov}+m$ is normal.  No code change; the earlier
  proposed DDH rewrite was reverted.

**Other.**

- Global scope of the study (multi-ensemble sweep, per-case current/mass handling) is now
  specified in Sec. 1.1 and Chunk 10.
- Tuning, not design: number of sink times $n_s$ and $\zeta$-noise hits for the diagonal
  projections (statistics vs cost); pin the $\kappa^{(0)}_{t,t+1}(n)$ accessor when coding.
- `apply_k_dag` is implemented and under verification; confirm before connected coding.
- Decide the output directory/filename scheme that encodes both sea and valence mass
  (Sec. 4) before coding Chunk 6/Chunk 10.
