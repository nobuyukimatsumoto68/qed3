# redo_obs_measurement_impl_plan_claude.md

Plan for measuring observables LOCALLY (barracuda22, 2x TITAN V + MPS) on the REDO production
ensembles generated at FNAL (affine, L1/L2 + L4 Nf4/6) and BU SCC (L4 Nf2). Ensemble inventory and
health: `redo_ensembles_claude.txt` (this directory; now also the obs-state blackboard, see Sec 6).

STATUS (2026-07-18): RUNNING. Configs local (rsync scripts); condensate DONE; conn ~87%;
glue shapes sweep running; disc NOT started; flow sweep at protocol v3 (Sec 4d) running.
Blackboard obs counts not yet filled. Session memory: project_redo_obs_flow.md.

## 1. Goal / observable matrix

| ensembles | obs | device | stride | priority |
|---|---|---|---|---|
| massless (9 L1 + 9 L2 FNAL, 3 L4 SCC Nf2, 6 L4 FNAL Nf4/6 not-launched) | conn jj+scalar $Y_{lm}$ tower (vector s1/s2/s3 + axial + $\sigma_{PS}/\sigma_{FS}$) | GPU | 20 | **1 (first GPU job)** |
| same massless | disc jj (vector+axial loops) | GPU | 20 | 2 |
| same massless | glueball corr: F shapes + F^2 shapes (7-shape consolidated basis) | CPU | 1 | parallel to GPU work |
| massive (4 L1 + 4 L2 FNAL, 4 L4 SCC) | condensate $\langle\sigma_{PS}\rangle,\langle\sigma_{FS}\rangle$ (e/o+spin dilution) | GPU | 20 (see Q5) | 3 (after conn) |
| ALL streams | NEW gluonic thermalization observable (Sec 4) | CPU | 1 | 0 (sets kmin for everything) |

- Fermionic stride = 20 (INCLUDING massive condensate, Q5 resolved), glue/therm stride = 1
  (NM 2026-07-17).
- Thermalization cut kmin: **SET k >= 20 (all massless, 2026-07-18).** From the flow monitor
  E_s(t0) at r/t0=1 vs trajectory (`therm_cut_first20_claude.png`): rise is k=1-11, plateau by
  k~12; k=20 = first config on the stride-20 jj grid, so the cut is free. Applied at ANALYSIS
  time (measurements stay uncut). Massive: TBD (short streams; dense-condensate monitor
  `run_cond_therm_claude.sh` written, unlaunched).
- Scalar DISC on massless: expected to VANISH at m=0 (dropped in the gsq8 campaign,
  scalar_ylm_corr plan 2026-07-15) -> default = conn scalar only; disc = jj vector+axial only.
  Confirm (Q2).

## 2. Source ensembles / rsync (Chunk 0)

Pull `ckpoint_lat.*` ONLY (no rng/logs/binaries), per the rsync recipes already in
`redo_ensembles_claude.txt` (FNAL block, lines at bottom; SCC block mid-file).

- FNAL: `nmatsum@lq.fnal.gov:/lustre2/affine/redo/` (L1/L2 all + L4 Nf4/6 when launched)
- SCC : `nmatsum@scc1.bu.edu:/projectnb/qfe/nmatsum/qed3/src/production/` (L4 Nf2, 7 streams)
- LOCAL DEST (proposal, Q3): `/mnt/barracuda22/qed3/qed3/src/production/<ensemble-basename>/`
  -- mirrors both remotes' layout, keeps redo separate from the legacy both_3d ensembles.
- Disk: total at full targets ~5 GB vs 4.3 TB free on /mnt/barracuda22 -> NO cleanup needed.
- Handoff scripts WRITTEN (2026-07-17, user runs them independently, repeatable as streams fill):
  `rsync_fnal_claude.sh` + `rsync_scc_claude.sh` (single rsync invocation each = one auth prompt;
  -z on, no -c; tee to rsync_{fnal,scc}_claude.log; per-ensemble ckpoint counts printed at the
  end for the blackboard loc_ncfg column; ssh ControlMaster hint in the SCC header).

## 3. Drivers to copy from both_3d (Chunk 1)

Copy (cp, not move) from `../both_3d/` into this directory; `includes/` here is already a superset
of `both_3d/includes/` (verified 2026-07-17, includes icos_orbits/wilson_shapes/flow headers), so
only the top-level `.cu` files move:

| driver | obs | notes |
|---|---|---|
| `jj_local_ylm_scalar_conn_stoch_claude.cu` | conn vector+axial+scalar tower | has `--ens-dir`; per-L binary (compile-time N_REFINE); FRESH mode (not `--IsScalarOnly`) |
| `jj_local_ylm_disc_stoch_claude.cu` | disc jj | has `--ens-dir`; tb2 self-contraction fix in place |
| `jj_scalar_condensate_eo_stoch_claude.cu` | condensate (massive) | has `--ens-dir`; e/o+spin dilution, 4 patterns/hit |
| `glue2_msm_shapes_claude.cu` | F shapes corr | dir3 built INTERNALLY from gsq/Nf/nu0 -> does NOT know the redo `mRe...mIm..._hb*` dir names -> needs Sec 5 fix |
| `glue_f2_shapes_claude.cu` | F^2 shapes corr | same dir3 issue |

Run-script templates to adapt (copy + edit, same scheduling ideas):
- `run_conn_gsq8_massless_claude.sh` -- GPU rotated-worker scheduler, config-range split, MPS
  auto-start, per-L binaries, resumable ("complete"-gated). Adapt: ensemble list from
  `redo_ensembles_claude.txt`, STRIDE=20, `--kmin` from therm cut, dest dirs.
- `run_interacting_shapes_claude.sh` -- CPU 8-worker glue sweep (OMP_NUM_THREADS=1 each).
  NOTE: when Claude runs CPU jobs directly, cap at 4 cores total (user rule); size the handoff
  script for the user's preference instead (ask at Chunk 4).
- `run_condensate_eo_heavy_claude.sh` -- condensate invocation pattern (auto-detect k range).

## 4. NEW gluonic thermalization observable (Chunk 2)

New CPU driver (name proposal: `glue_therm_series_claude.cu`), one row per config -> MC-time
series used to set the per-ensemble thermalization cut kmin.

SCHEME (NM 2026-07-17, superseding the fixed-t_fl proposal): use the WILSON-FLOW SCALE $t_0$,

$$t_0^2 \, E_s(t_0) = c , \qquad c = 0.3 \ \text{(first guess; tuned after the trajectory check)}$$

- $E_s(t_{fl})$ = SPATIAL average action density, SAME $\beta_s$ weight as the glueball code
  (noncompact $0.5\beta_s\theta^2$ integrand per face, l=m=0 of density_Ylm_spatial; no temporal
  channel in the observable).
- FLOW = per-timeslice SPATIAL gradient flow of the simulated action (get_spatial, $\beta_s$
  weight; NOT coupling-stripped). SWITCHED from full 3D (NM 2026-07-17 later the same day):
  matches the glueball "on a timeslice" convention and avoids the anisotropic two-rate kernel
  ($\beta_s \sim a_t$ vs $\beta_t \sim 1/a_t$) -- spatial-only has ONE uniform prefactor
  $a_t/g^2$, so $a_t$ drops out of the $t/g^2$ collapse bookkeeping. LO scaling unchanged
  ($E_s \sim g^2/t^{3/2}$; the $p_0$ integral converges without temporal damping).
- Driver records the TRAJECTORY $E_s(t_{fl})$ every SAVE_EVERY=20 steps: $t_{fl}=0,0.2,...,2.0$
  (dt=0.01, tmax=2.0; first-try values, tunable per ensemble). Smooth curve + $t_0$
  reconstructed in analysis; the MC series of $E_s$ at the chosen $t_0$ is the therm monitor.
### 4b. Renormalized coupling + scale from the flow (NM question 2026-07-17; CORRECTED)

LO with the $[B_\mu]=1$ convention ($[E]=4$ in any $d$; $[g^2]=1$ in 3D): the flowed
propagator carries $e^{-2tp^2}$, so

$$
\langle E(t)\rangle_{LO} \;=\; C\, g_0^2 \int \frac{d^3p}{(2\pi)^3}\, e^{-2tp^2}
\;\propto\; \frac{g_0^2}{t^{3/2}} .
$$

NOTATION (NM question 2026-07-17): in the LO formula above, $g^2$ = the BARE coupling
$g_0^2$ (the ensemble input gsq) -- at LO bare vs renormalized is indistinguishable. The
renormalized coupling $g_R^2 \equiv g^2_{GF}(t)$ is DEFINED by reading the relation
backwards, $g^2_{GF}(t) = t^{3/2}\langle E_s(t)\rangle/C'$: the flow renders $E$ finite, so
this is a legitimate nonperturbative definition (scheme: GF, spatial-$E$, per-timeslice
flow). In superrenormalizable 3D the coupling needs no infinite renormalization; the finite
corrections organize in the dimensionless parameter $g_0^2\sqrt{t}$:
$g^2_{GF}(t) = g_0^2(1 + k\,g_0^2\sqrt{8t} + \dots)$. Hence (i) the plateau is only
LO-flat -- the $\sqrt t$ NLO growth is visible as the plateau-height/$g^2$ increase with
coupling in the v3 check data; (ii) $R$ = the $S^2$ sphere radius ($=1$), the IR scale
bounding the flow window from above (flow radius $\ll R$), as $a$ bounds it from below.

ADOPTED DEFINITIONS (NM 2026-07-17, v6 protocol): $t_0$ from $t_0^2 E_s(t_0) = c$ with
$c = 0.003$ (updated from 0.002, NM 2026-07-17); renormalized coupling $g_R^2 \equiv t_0^{3/2}E_s(t_0)$. These are ONE
measurement: the identity $g_R^2 = c/\sqrt{t_0}$ holds exactly. LO guideline:
$t^{3/2}E_s = C' g_0^2$ (flat) $\Rightarrow t_0^{LO} = (c/(C'g_0^2))^2$,
$g_R^{2,LO} = C'g_0^2$; the natural x-axis of the $t^{3/2}E_s$ plot is $t/t_0^{LO}
\propto t\,g_0^4$ (all LO crossings at $x=1$; horizontal shift of the real crossing =
$g_R^2/g_0^2$). $C'$ from the weakest-coupling flat regions (r=0.5 trio ~ free limit).

### 4c. Lattice-unit bookkeeping: what $\hat E_s$ actually is (RESOLVES the "$\times L$" question, 2026-07-17)

The driver stores the PURE LATTICE NUMBER
$\hat E_s = \langle \tfrac12\beta_s\theta^2\rangle_{faces}$ with $\beta_s = a_t/(g_0^2 V_i)$,
$\theta = V_i F_{12}$. Per face this is $a_t V_i\, E_B(x)$ with $E_B \equiv F_{12}^2/(2g_0^2)$,
the ACTION-DENSITY convention ($[E_B]=3$, coupling INSIDE -- not the $[E_A]=4$ convention of
4b, $E_A = g_0^2 E_B$). Mean face area on the unit sphere: $\bar V = 4\pi/(20L^2) = \pi/(5L^2)$.
Hence with the flow-time conversion $t \propto a_t\hat t/g_0^2$:

$$
t^{3/2} E_B \;\propto\; \hat t^{3/2}\hat E_s \cdot \frac{\sqrt{a_t}\, L^2}{g_0^3} ,
\qquad
t^{3/2} E_A \;\propto\; \hat t^{3/2}\hat E_s \cdot \frac{\sqrt{a_t}\, L^2}{g_0} .
$$

At LO $t^{3/2}E_B$ = pure diffusion constant (NO coupling) $\Rightarrow$ the LAB plateau scales
as $\hat t^{3/2}\hat E_s \propto g_0^3/L^2$. DATA CHECK (v6 plateaus): L2 couplings 1:2:3 give
5.54e-3 : 1.59e-2 : 2.98e-2 = 1 : 2.87 : 5.38 vs $g_0^3$ prediction 1 : 2.83 : 5.20 (the naive
$g_0^2$ would be 1:2:3 -- clearly excluded); fixed $g_0^2=2$: L2/L4 plateau ratio 4.7 vs $L^2$
ratio 4. CONFIRMED. Consequences:
- the "$\times L$" hunch was pointing at this conversion; the exact factor is $L^2/g_0$ for
  the COUPLING readout: $g_R^2 \propto \hat t^{3/2}\hat E_s(\hat t_0)\, L^2/g_0$ (reduces to
  $\propto g_0^2$ at LO), with one overall constant calibrated so $g_R^2 \to g_0^2$ at weak
  coupling;
- the UNIVERSAL (cross-ensemble) dimensionless scale condition is
  $\hat t^2 \hat E_s \cdot L^2/g_0^2 = \hat c$ (equals const $\times\, t_{phys}^2 E_A$);
  the current per-L $C'_L$ device was implicitly absorbing exactly this $g_0/L^2$ factor
  (measured $C'_L$ ratios 1 : 0.46 : 0.15 vs the $g_0/L^2$ prediction 1 : 0.35 : 0.13,
  residual = renormalization).

TWO distinct objects follow (the first 4b draft conflated them):

1. **Renormalized coupling (mass units)**:
   $t^{3/2}\langle E(t)\rangle = C'\, g_0^2 + O(g_0^4)$ has mass dimension 1 and is
   $t$-INDEPENDENT at LO. It is the GF-scheme $g^2_R$ read off as a PLATEAU in the scaling
   window $a^2 \ll t_{phys} \ll R^2$; the small-$t$ rise of our curves is the discretization
   region (flow radius ~ lattice spacing), the large-$t$ decay is the finite sphere. A flat
   region between them = the measurement of $g^2_R$. (Window is narrow/absent at L1,
   opens with L -- visible in the check data: L4 curves are flattest near their peak.)
2. **Scale setting (dimensionless condition)**: the dimensionless running coupling at
   $\mu = 1/\sqrt{8t}$ is
   $$
   \hat g^2(\mu) = \frac{g^2_{GF}(t)}{\mu} = \frac{\sqrt{8t}}{C'}\; t^{3/2}\langle E\rangle
   \;\propto\; t^2 \langle E(t)\rangle ,
   $$
   i.e. the 4D-LOOKING combination $t^2\langle E\rangle$ is the correct dimensionless one in
   3D too -- it grows $\propto g_0^2\sqrt{t}$ CLASSICALLY (superrenormalizable: powerlike, not
   log, running). The condition $t_0^2\langle E(t_0)\rangle = c$ then gives
   $\sqrt{8 t_0}\, g_R^2 \approx c/C'$, so $1/\sqrt{t_0} \propto g_R^2$: $t_0$ tracks the
   physical correlation scale set by the renormalized coupling, and the crossing is
   WELL-CONDITIONED (linear growth in $\sqrt t$).

CONCLUSION: use $t^2 E_s = c$ for the scale $t_0$ (c tuned to the measured range, ~few
$\times 10^{-3}$ from the Nf2 check); use the $t^{3/2}E_s$ plateau for $g^2_R$. A
$t^{3/2}E_s = c$ crossing on its steep rising branch would sit in the artifact region
(flow radius ~ a) and would NOT define a continuum-legit scale -- this was the flaw in the
first draft. Same three caveats as below still apply:

1a. **Lattice -> continuum flow-time translation (NM derivation 2026-07-17, D-corrected)**:
   continuum $S = \frac{1}{g^2}\int d^Dx\, F^2$, $[A]=1$, $[g^2]=4-D=1$ (D=3);
   $U \sim 1 + aA$. Lattice Wilson flow (everything dimensionless):

   $$\frac{dU}{d\hat t}U^{-1} = -\partial_U S
   \;\Longrightarrow\; a\,\frac{dA}{d\hat t} = -a^{D-1}\,\frac{\delta S}{\delta A}
   \;\Longrightarrow\; \frac{dA}{d\hat t} = -a^{D-2}\,\frac{\delta S}{\delta A},$$

   using $\partial_\theta S = a^{D-1}\,\delta S/\delta A$ ($\theta = aA$, one site volume
   $a^D$ per link). The identification "$t=a^2\hat t$, $dA/dt=-\delta S/\delta A$" closes
   ONLY in D=4 (where $a^{D-2}=a^2$ and $[g^2]=0$). In D=3 the same manipulation leaves
   $dA/d\hat t = -a\,\delta S/\delta A$ -- dimensionally forced, since
   $[\delta S/\delta A]=D-1=2$ while a flow in $t=a^2\hat t$ would need dimension 3.
   Matching to the canonical Luscher flow $dA/d\tau = D\!\cdot\!G$ (no $1/g^2$;
   $\delta S/\delta A = \frac{2}{g^2}D\!\cdot\!G$) gives the physical flow time

   $$\tau \;=\; \frac{2\,a^{D-2}}{g^2}\,\hat t \qquad (D=3:\ \tau = 2a\,\hat t/g^2).$$

   ANISOTROPIC case (ours): site volume $a_s^2 a_t$, spatial-link variation removes one
   $a_s$: $\partial_\theta S = a_s a_t\,\delta S/\delta A$, so $dA/d\hat t = -a_t\,\delta
   S/\delta A$ and $\tau \simeq 2 a_t \hat t/g^2$ -- the prefactor is $a_t$ (NOT $a_s$),
   independently confirming the set_beta result below including which lattice spacing
   appears. Since $a_t=0.2$ is COMMON to all ensembles, the cross-ensemble collapse
   variable remains $\hat t/g^2$.

1b. **Flow-time normalization (DERIVED 2026-07-17 from set_beta, action_ext_claude.h)**:
   $\beta_s[i] = a_t/(g^2 V_i)$ with face area $V_i\sim a^2$, and
   $\beta_t \sim 2/(a_t g^2)$ ($V_{link}/\ell^2 = O(1)$). Linearizing the noncompact flow,
   a physical-momentum-$p$ mode damps at rate $\sim (a_t/g^2)\,p^2$ per unit lab flow time:

   $$t_{phys} \;\simeq\; \frac{a_t\, t}{g^2} \quad (\text{up to } O(1) \text{ geometry; L-independent}).$$

   COLLAPSE VARIABLE = $t/g^2$, NOT the guessed $t/r$ (no $L$ factor). Coarse-grid peaks
   support this: $t_{peak}/g^2 \approx$ const within each $L$ (0.2 / 0.07 / 0.02 for
   L1/L2/L4) and that constant scales like $a^2 \sim 1/L^2$ -> the PEAK is the flow radius
   reaching ~one lattice spacing (artifact boundary), the tail beyond is the continuum-like
   window (widest at L4). Prediction for the fine-grid check: same-$L$ curves stack in
   $t/g^2$; cross-$L$ tails overlap; peaks sit at $t/g^2 \propto 1/L^2$. Analysis plots BOTH
   $x = t/r$ (NM request) and $x = t/g^2$ (derived) to let the data arbitrate.
2. **Finite sphere**: $R=1$ is an IR scale; the peak/turnover of $t^{3/2}E_s$ is the flow
   radius saturating the sphere. Defining the coupling on the RISING branch (radius << R)
   approximates the infinite-volume scheme; alternatively embrace a finite-volume scheme with
   $\sqrt{8t_{phys}}/R$ fixed (4D GF-torus style).
3. **Normalization $\mathcal{N}$ / scheme spec**: $E_s$ = spatial (magnetic) term only, on an
   anisotropic curved lattice -- a valid scheme, but $\mathcal{N}$ differs from the full-$E$
   LO constant. Get $\mathcal{N}$ nonperturbatively from the FREE-FIELD flow (same driver,
   free config): there $t^{3/2}E_s/g_0^2$ at small coupling gives $C'$ directly, so
   $\hat g^2_{GF}$ matches $g_0^2$ at weak coupling.

IMPLEMENTED: `glue_therm_series_claude.cu` (this dir; SPATIAL flow (no FLOW_FULL), explicit
ens_dir arg 7, output `data_<basename>/therm_series_claude.dat` TEXT append+skip resumable;
production h5 layout = PLANNED LATER). Handoff `tmp_glue_therm_claude.sh`: build L1/L2/L4
binaries (-DN_REFINE_CLI) + FIRST TRAJECTORY CHECK on 6 local redo ensembles (Nf2 weak/strong
gsq per L; L1 stride 100, L2 stride 50, L4 stride 2) to set c. Cut rule: kmin = first k where
the $E_s(t_0)$ series is inside ~1 sigma of the tail plateau (visual check per ensemble);
record kmin per ensemble in the blackboard and use as `--kmin` for the fermionic obs.

### 4e. Curved-kernel LO (S^2 x R mode sum) and the sqrt(L) collapse (2026-07-18)

On $S^2\times R$ the flat $\int d^3p$ of Eq. (1) is replaced by the transverse gauge-mode sum
of qed3_v2-6.pdf App. C (Eq. C.37): eigenvalues $\lambda = k^2 + j(j+1)$, $j = n+|m| \ge 1$
(curvature gives the transverse modes a "mass" $j(j+1)/R^2$), degeneracy $\propto (2j+1)$.
The per-timeslice spatial flow damps only the spatial part: mode $j$ decays at rate
$\kappa\,\lambda_j$ per unit lab time, $\lambda_j = j(j+1)$, $\kappa \propto a_t/g_0^2$.
Doing the (undamped) $k$ integral, $\int \frac{dk}{2\pi}\frac{\lambda_j}{k^2+\lambda_j} =
\frac{\sqrt{\lambda_j}}{2}$:

$$
E_B^{LO}(\hat t) \;\propto\; \sum_{j\ge1} (2j+1)\,\sqrt{j(j+1)}\;
e^{-2\kappa\, j(j+1)\,\hat t}\,, \qquad \kappa \propto \frac{a_t}{g_0^2}\,.
$$

Regimes:
1. WINDOW ($a \ll$ flow radius $\ll R$): sum $\to \int dj\, j^2 e^{-2\kappa j^2\hat t}
   \propto (\kappa\hat t)^{-3/2}$ -- the flat-space law. In lattice variables
   ($\hat E_s \propto E_B/L^2$): $\hat E_s \propto g_0^3/(L^2\hat t^{3/2})$, and AT FIXED
   $r = g_0^2/L$ (so $g_0^3 = (rL)^{3/2}$):
   $$\hat E_s\sqrt{L} \;=\; C\,(r/\hat t)^{3/2} \quad\text{-- the observed sqrt(L) collapse}$$
   is the CONTINUUM-WINDOW PREDICTION at fixed $r$, NOT a coincidence. Checks (h5 data,
   $\hat t = 2$): family ratios 1.97 / 2.75 vs predicted $(1.5)^{3/2} = 1.84$ /
   $2^{3/2} = 2.83$; decay $\hat t\,0.5\to2$: factor 7.5 vs 8.
2. CURVATURE/IR ($\kappa\hat t \gtrsim 1$): sum truncates to $j=1$:
   $\hat E_s \to$ const $\times\, e^{-4\kappa\hat t}$, exponential with rate
   $4a_t/g_0^2 \propto 1/(rL)$ at fixed $r$ -- the collapse BREAKS $L$-dependently
   (weakest/finest first), as the data show.
3. CUTOFF (small $\hat t$): $j_{max} \sim$ few$\times L$ truncation = artifact region.

### 4f. Why $g_0^2\bar a_s \sim g_0^2 R/L$ is the perturbative coupling (NM question 2026-07-18)

Conventions from qed3_v2-6.pdf: $R=1$ units, spatial spacing $\bar a_s \propto R/L$ ("$\bar a_s$
defines the cutoff scale for the state profile on the sphere", Sec. II); gauge action
$\frac{1}{4g_0^2}\int F^2$ (lattice: $\beta_s = a_t/(g_0^2 V_i)$), transverse mode eigenvalues
$\lambda = k^2 + j(j+1)$, $j\ge1$ (Eq. C.37).

1. With this normalization the PHOTON propagator carries the coupling,
   $\langle AA\rangle \sim g_0^2/\lambda$, while vertices and the fermion propagator are
   $g_0^2$-independent. Adding one internal photon line to any diagram therefore multiplies
   it by $g_0^2$ -- and since $[g_0^2]=1$ and nothing else in the loop is dimensionful, by
   one INVERSE power of the characteristic virtuality $p$ of that line:

   $$\text{extra photon line} \;\Rightarrow\; \times\; \frac{g_0^2}{p}\,.$$

2. MODE-RESOLVED coupling on the sphere ($p \sim j/R$):
   $$\lambda_{eff}(j) \;=\; \frac{g_0^2 R}{j}\,, \qquad
   \lambda_{eff}(1) = g_0^2R \ (\text{IR strength}), \quad
   \lambda_{eff}(j_{max}\sim L) \sim g_0^2\bar a_s \sim \frac{g_0^2R}{L} \sim r\,.$$
   ONE theory contains BOTH regimes: the sphere-scale modes are strongly coupled when
   $g_0^2R \gg 1$, while the cutoff modes are coupled with strength $r$.

3. UV divergences of any diagram come exclusively from the $p \sim 1/\bar a_s$ region, where
   each additional order costs $g_0^2\bar a_s$. Hence higher orders are MORE UV-convergent
   (the superficial divergence degree DECREASES with the order) -- this is
   superrenormalizability restated -- and the counterterm structure is an expansion in
   $g_0^2\bar a_s \to 0$, EXACT in the continuum limit at ANY fixed $g_0^2R$. Strong IR
   coupling never enters the UV analysis (locality/decoupling).

4. Concrete check with the paper's objects: the one-loop vacuum polarization in 3D is
   $\Pi(p) \simeq N_f\,p/16$ ($g_0$-independent in this normalization); it corrects the
   inverse propagator $\lambda/g_0^2 + \Pi$, i.e. a RELATIVE correction
   $g_0^2\Pi/\lambda \sim N_f g_0^2/(16 p)$: at $p\sim1/\bar a_s$ this is
   $N_f g_0^2\bar a_s/16$ (small as $a\to0$); at $p \sim 1/R$ it is $N_f g_0^2 R/16$
   (the strong IR screening). Same diagram, same parameter $g_0^2/p$, two regimes.

CONSEQUENCE for our trajectories: $r = g_0^2/L \in [0.5, 1.5]$ means the CUTOFF modes are
currently O(1)-coupled -- cutoff corrections at fixed $r$ are not perturbatively small and
the continuum limit must be taken at FIXED $g_0^2$ (r decreasing), where they vanish.
Data-level nonperturbative support: the $E_sL^2$ matching of equal-$g_0^2$ pairs (L2/L4 g2,
L1/L2 g1) with NO recalibration of $g_0^2$.

### 4d. PRODUCTION flow driver (WRITTEN 2026-07-17, protocol frozen)

`glue_therm_flow_claude.cu`: per-ensemble HDF5 `data_<basename>/therm_flow_claude.h5` with
`tlist` [201] ($\hat t = 0, 0.01, \dots, 2.0$; eps 0.01, save every step, tmax = arg 8),
`E` [n_t][n_cfg] (config = fast axis, jackknife per fixed $\hat t$), `klist` sorted, and grid
metadata. RESUME reads the file and flows ONLY missing configs (flow never recomputed);
rewrite after every new config via tmp+atomic-rename; grid-mismatch ABORTS (no mixing).
Build+smoke handoff: `tmp_glue_therm_flow_claude.sh` (3 binaries, L1 g0.5 stride-100 smoke,
resume-verification rerun, h5 shape dump). PENDING: all-ensemble sweep script (stride 1,
massless + massive, 8-worker pool); therm-cut/t0/g_R analysis on the h5 (tomorrow's
normalization discussion first -- Sec 4c).

## 5. Glue drivers vs redo dir names (Chunk 3)

Redo dirs are `Nf*_gsq*at*nu0*mRe*mIm*nt128L*_hb*` -- glue drivers' internal dir3 construction
(`glue2_msm_shapes_claude.cu:246-256`) cannot form these. DECIDED (Q4, NM 2026-07-17): PATCH the
copied drivers with an explicit ensemble-dir argument (original dir3 line left commented per
migration convention), and write the output h5 into the SAME `data_<ensemble-basename>/` dir the
fermionic obs use (fermionic drivers write to `data_<ESNID>/...`).
DONE 2026-07-17: both copies patched (arg7 = ens_dir, data_<basename>/ output). Sweep handoff
`run_glue_shapes_redo_claude.sh`: all local MASSLESS ensembles, EVERY config (stride 1, kmin 1),
7-shape basis with the phase-2 saved l-towers KEPT (F l=1,2; F^2 l=0,1 -- NM: no drop-op pruning,
towers as saved), 8 single-threaded workers (xargs pool, largest ensemble first, NWORK override),
resumable via complete-gated h5. Runs alongside the GPU campaign `run_cond_conn_claude.sh`
(4 workers, 2/GPU under MPS -- NM corrected: MPS packing IS used).

## 6. Blackboard obs-state columns (DONE in this prep pass)

`redo_ensembles_claude.txt` rows gain two trailing columns (mirroring both_3d
`jj_ensembles_claude.txt` `existing_jj` style):
- `loc_ncfg` -- ckpoint_lat configs pulled to barracuda22 so far (0 = not pulled yet);
- `obs` -- comma-joined `tag:N` (N = configs measured locally), `-` = none. Tags:
  `thm` (therm series) | `conn` | `disc` | `glue` (F) | `f2` (F^2) | `cond` (condensate).
Update these counts as pulls/runs progress (same discipline as the both_3d manifest refresh).

## 7. Open questions (NM) -- Q1-Q5 RESOLVED 2026-07-17

- **Q1 therm observable** RESOLVED: plaquette + flowed plaquette, FULL 3D flow (Sec 4).
- **Q2 scalar disc** RESOLVED: stays dropped; massless disc = jj vector+axial only.
- **Q3 local dest** RESOLVED: `src/production/<basename>/` (rsync scripts written, Sec 2).
- **Q4 glue dir fix** RESOLVED: patch drivers; output into the fermionic `data_*` dirs (Sec 5).
- **Q5 massive stride** RESOLVED: 20 (same as all fermionic obs).
- **Q6 conn/disc params** OPEN: keep the gsq8-campaign settings otherwise (nhits=1, t0=0,
  spin-dilution ON, l=0..3, tb2 for disc)? Condensate: valence mass = sea mRe, nhits as in the
  heavy runs (nhits=1)?
- **Q7 L4 Nf4/6** OPEN: FNAL slots p13-p15 are not-launched -- launch is a separate (remote)
  action; local measurement just picks them up via rsync once configs appear. Confirm nothing to
  do locally now.

## 8. Chunk order

0. Rsync handoff script + first pull (needs Q3).
1. Copy drivers + build all per-L binaries (build handoff or user-run make; no runs yet).
2. Therm driver: write + run stride 1 on everything pulled (needs Q1) -> set kmin per ensemble.
3. Glue driver dir fix (needs Q4) + glue F/F^2 CPU sweep stride 1 (runs alongside GPU jobs).
4. Conn GPU campaign stride 20 (priority 1), then disc (2), then massive condensate (3, needs Q5).
5. Blackboard count refresh after each stage.
