# Plan v2: Conserved-Current Correlator Measurements

Clean re-summary of the design worked out in
`conserved_current_correlators_impl_plan_v1_claude.md` (v1, scratch).  All decisions below
are settled; v1 keeps the derivation history.

Reference: A. Katz, N. Matsumoto, "Notes on QED3 (interacting theory)",
`/mnt/barracuda22/qed3/qed3int_v2.pdf` -- Secs. 3.3 (lattice correlators, Eqs. 3.39--3.57),
4.2--4.3 (CFT estimators, Eqs. 4.27--4.36), App. C (real spherical harmonics).

---

## 1. Goal and scope

Measure current-current two-point functions on the $S^2\times$time lattice using the
exactly conserved overlap currents, for $N_f = 2, 4, 6$.  Currents are two-component in the
$\xi$/$\eta$ basis; the (+) vector current is $J_{V,+}^{nn'}(t) = \eta_f^\dagger K^{nn'}(t)\xi_f$.

Current types:

- **Vector** $J_{V,\pm}$ -- connected $+$ disconnected (Eqs. 3.41--3.44).
- **Flavor** $J_{r,\pm}$ -- connected only, equal to the vector-connected piece
  (Eqs. 3.49--3.52); **does not exist for $N_f = 2$** (relevant for $N_f = 4, 6$).
- **Axial** $J_{A,\pm}$ -- connected only, both orderings $C_A(0\to t)$, $C_A(t\to 0)$
  (Eqs. 3.45--3.48).

Three projected observables, computed for each current type:

- $G_s(t)$ -- spatial, equal-link average (Eq. 4.29).
- $G_t(t)$ -- temporal, site-diagonal (Eq. 4.32).
- $G^{t}_{\ell_1 m_1;\ell_2 m_2}(t)$ -- spherical-harmonic projection (Eq. 4.36), which
  carries the conformal tower (Eqs. 4.34--4.35): the $\ell=0$ channel vanishes (charge
  conservation), $\ell=1 \sim e^{-2t}$, $\ell=2 \sim e^{-3t}$.

---

## 2. Notation

Fermionic contractions (Eqs. 3.39--3.40), with the fermionic average $\langle\cdot\rangle_F$:

$$
\langle \xi_{f,x_1}\,\eta_{g,x_2}^\dagger\rangle_F = \delta_{fg}(D_\text{ov}^{-1})_{x_1 x_2},
\qquad
\langle \eta_{f,x_1}\,\xi_{g,x_2}^\dagger\rangle_F = \delta_{fg}(D_\text{ov}^{-\dagger})_{x_1 x_2}.
$$

There is **no $\gamma_5$-hermiticity** here ($D_\text{ov}^{-\dagger}\neq\gamma_5 D_\text{ov}^{-1}\gamma_5$).
$\langle\cdot\rangle_U$ is the gauge average; $^*$ is complex conjugation.

Kernels: $K^{nn'}(t)$ (spatial link) and $K^{t,t+1}(n)$ (temporal link), via
`ConservedCurrent::apply_k`; the adjoint $K^\dagger$ via `apply_k_dag`
(`conserved_current_claude.h:157`, spatial wrapper `apply_k_dag_wz` `:230`; implemented,
under verification).  Both accept `std::pair<int,BaseLink>` (spatial) and
`std::pair<int,Idx>` (temporal).

---

## 3. Estimators

### 3.1 Disconnected (vector only) -- factorizes

$$
C_{V,d}^{nn'}(0\to t) = T^{nn'}(0)\,T^{nn'}(t),
\qquad
T^{nn'}(t) \equiv \mathrm{tr}[D_\text{ov}^{-1}K^{nn'}(t)] \approx \eta^\dagger K^{nn'}(t)\,\phi,
\quad \phi = D_\text{ov}^{-1}\eta,
$$

estimated with $Z_2\times Z_2$ noise and time-spin dilution (as `disc_claude.cu`).  The
measurement program **dumps the raw per-config single-time traces** $T^{nn'}(t)$ (spatial)
and $T^{t,t+1}(n)$ (temporal); the product, gauge average, vacuum subtraction, weighting,
and projection are done later (Secs. 3.4--3.5, 5).

### 3.2 Connected (vector / flavor) -- single source, source/sink split

$$
\phi = K^{(\text{src})}\,D_\text{ov}^{-1}\eta,
\qquad
\psi = D_\text{ov}^{-\dagger}\,K^{(\text{snk})\dagger}\,\eta,
\qquad
C_{V,c} = \langle \psi^\dagger\phi\rangle .
$$

Since $\psi^\dagger = \eta^\dagger K^{(\text{snk})}D_\text{ov}^{-1}$,
$E[\psi^\dagger\phi] = \mathrm{tr}[D_\text{ov}^{-1}K^{(\text{snk})}D_\text{ov}^{-1}K^{(\text{src})}] = C_{V,c}$:
a single bilinear, exact, both legs forward $D_\text{ov}^{-1}$, no disconnected contamination.

Efficiency: the **source side is free** -- $\phi_0 = D_\text{ov}^{-1}\eta$ is solved once
per hit and the source kernel is applied after it, covering all source insertion points and
times.  The **sink side** is one $D_\text{ov}^\dagger$ solve per sink time on a *composite*
source (Sec. 3.4).

### 3.3 Axial (both orderings) -- same split with GW factors

Forward (Eq. 3.45),
$C_A(0\to t)=\mathrm{tr}[D_\text{ov}^{-1}K(0)(1-D_\text{ov}^\dagger)D_\text{ov}^{-\dagger}K^\dagger(t)(1-D_\text{ov})]$:

$$
\phi_A = D_\text{ov}^{-\dagger}\,K^\dagger(t)\,(1-D_\text{ov})\,\eta \ (\text{sink side}),
\qquad
\psi_A = (1-D_\text{ov})\,K^\dagger(0)\,\rho_0 \ (\text{source side, free}),
\qquad \rho_0 = D_\text{ov}^{-\dagger}\eta,
$$

with $C_A(0\to t) = \langle\psi_A^\dagger\phi_A\rangle$.  Backward $C_A(t\to 0)$ (Eq. 3.47)
is the $K\leftrightarrow K^\dagger$, $D_\text{ov}\leftrightarrow D_\text{ov}^\dagger$ mirror
(source side uses $\phi_0=D_\text{ov}^{-1}\eta$).  Implemented **as written**.

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

For the vector current $\langle JJ\rangle = -C_{V,c}+C_{V,d}$.  The geometric weights (and
the real $Y_{\ell m}$ from `GetYlm`) are applied **inside the compiled programs**; the
notebooks only gauge-average, vacuum-subtract, combine, and fit.

**Normalization convention.**  The $C$'s above are the bare traces; the PDF relations carry
overall factors, e.g. $\tfrac{2}{N_f}\langle J_{V,+}(0)J_{V,+}(t)\rangle_F = -C_{V,c}+C_{V,d}$
(Eq. 3.41) and the disconnected loop is $L=-N_f\,T$.  These overall $2/N_f$, $-N_f$ factors
are bookkeeping applied at the analysis stage; the programs output the bare $C_{V,c}$,
$C_A$, $T$.

### 3.6 Geometric quantities (already in code)

| Symbol | Meaning | Accessor |
|---|---|---|
| $\kappa^{(0)}_{nn'}$ | free-field spatial-link coupling | `DiracOp::kappa[il]` (`dirac_simp.h:355-360`) |
| $A_{nn'}$ | dual area of edge $(n,n')$ | `lattice.link_volume[il]` (`s2n_simp.h:29`) |
| $A_n$ | dual (Voronoi) area of site $n$ | `lattice.dual_areas[s]` (`s2n_simp.h:38`) |
| $\kappa^{(0)}_{t,t+1}(n)$ | temporal-link coupling (uniform in $t$, per site) | per-site array (candidate `dual_areas[ix]`); pin when coding |
| $Y_{\ell m}(\hat n)$ | real spherical harmonic | `lattice.GetYlm(s,l,m)` after `UpdateYlm(l_max)` (`s2.h:55-56`) |

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
- Production noise: `nhits=1`, `t_block=8` (time-spin dilution), as `disc_claude.cu:164-165`.
  Connected programs use all-to-all $\eta$ (no `--t_block`); diagonal projections add
  $\zeta$ insertion-point noise.
- Output: one HDF5 file per config $k$, `HighFive::File::Truncate`, with resume sentinel
  (as `disc_claude.cu:272-286`).

Compile-time (inherit `disc_claude.cu`): `N_REFINE=1`, `Nt=128`, `NPARALLEL_DUPDATE=1`,
`TOL_INNER=1e-15`, `TOL_OUTER=1e-14`.

Reused infrastructure (not modified): `apply_k` / `apply_k_dag`
(`conserved_current_claude.h`); CG solve chain `op_DH` / `op_DHD` and a $D_\text{ov}^\dagger$
solve `op_DdagInv` (analogous normal-equation chain with $D_\text{ov}D_\text{ov}^\dagger$),
`from_cpu` / `solve` (`disc_claude.cu:231-237,294-295`); `FermionVector` with
`fill_z2_source`, `time_spin_dilution`, `dag` (`valence_claude.h`); `GaugeExt::read`
(`disc_claude.cu:282`); lattice geometry (Sec. 3.6).

---

## 5. Implementation chunks

### Disconnected dump

**Chunk 1 -- boilerplate + construction.**  Copy `disc_claude.cu` boilerplate (includes,
`Comp`, type aliases, `ParseArgs`, GPU/lattice/DW/Overlap/Gauge setup); add
`#include "includes/conserved_current_claude.h"` and `ConservedCurrent<Fermion> kop(D)`.
Files: `jj_disc_claude.cu`.

**Chunk 2 -- enumerate + dump layout.**  Enumerate spatial links $(s, lk)$ over `base.links`
and temporal sites $(n, s)$; HDF5 keys `spatial/{ix0}/{ix1}/{real,imag}` (len $N_t$) and
`temporal/{n}/{real,imag}` (len $N_t$); output dir + resume sentinel.
Files: `jj_disc_claude.cu`.

**Chunk 3 -- stochastic loop + dump.**  Per config, per hit, per dilution block: draw
diluted $\eta$, solve $\phi=D_\text{ov}^{-1}\eta$, accumulate
$T^{lk}(s)\mathrel{+}=\eta^\dagger K^{lk}(s)\phi$ and $T^{n}(s)\mathrel{+}=\eta^\dagger K^{t,t+1}(n,s)\phi$;
divide by `nhits`; write raw $T$ to HDF5.  No correlator formed.
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
(vector: $C_{V,c}$; axial: $C_A(0\to t)$ and $C_A(t\to 0)$).  Average over $\eta$ (and
$\zeta$) hits.
Files: the three `jj_conn_*proj_claude.cu`.

**Chunk 6 -- output.**  HDF5 per config; filename encodes `--current`; datasets: sp/tp ->
$G(t)$ length $N_t$; ylm -> per $(ch_1,ch_2)$ length $N_t$; axial writes two ($0\to t$,
$t\to 0$).  Resume sentinel.
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

---

## 6. Pending / future

- **Further complications** to the study -- to be specified (user); no rush.
- Tuning, not design: number of sink times $n_s$ and $\zeta$-noise hits for the diagonal
  projections (statistics vs cost); pin the $\kappa^{(0)}_{t,t+1}(n)$ accessor when coding.
- `apply_k_dag` is implemented and under verification; confirm before connected coding.
