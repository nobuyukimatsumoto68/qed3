# Deriving the overlap-fermion Hamiltonian and the energy density $T_{00}$

Goal: define $T_{00}$ PRECISELY as the conserved energy density = Hamiltonian density, from the lattice
action and exact temporal-lattice translation invariance -- replacing the naive
$\bar\psi\sigma_3\partial_0^{\rm sym}\psi$ interpolator, which we now know is only an approximation.

## Action (to confirm)
$$
S = \eta^H D_{\rm ov}\,\xi + \xi^H D_{\rm ov}^H\,\eta = 2\,\mathrm{Re}\big(\eta^H D_{\rm ov}\,\xi\big),
$$
Hermitian. Open: are $\eta,\xi$ independent Grassmann fields (the signed-mass pair, $\xi\!:\!+m$,
$\eta\!:\!-m$), and is $D_{\rm ov}$ the $S^2\times R$ overlap with $\sigma_3$ as the temporal gamma?

## Principle: the energy density generates the transfer matrix
- Exact symmetry: $S$ is invariant under the unit temporal shift $t\to t+1$ (fields AND links shifted).
  Noether $\Rightarrow$ a conserved energy.
- Transfer matrix $\hat T$ implements one time step: $Z = \mathrm{Tr}\,\hat T^{N_t}$, and
  $$
  \hat H = -\tfrac{1}{a_t}\log \hat T,\qquad \hat T\,|n\rangle = e^{-a_t E_n}\,|n\rangle .
  $$
- Any interpolator's correlator is $C_O(t)=\sum_n |\langle 0|O|n\rangle|^2\, e^{-a_t E_n t}$. The STATE
  ENERGIES $E_n$ are exact $\hat T$ data, independent of $O$. So for state ID any $0^{++},\ell=0,\Delta=3$
  operator suffices; the PRECISE $T_{00}$ is needed only for the physical density (Ward identity,
  normalization, protected dimension), not for reading $E_n$.

## The overlap obstruction (why the naive operator is imprecise)
$D_{\rm ov} = 1 + \gamma_5\,\mathrm{sign}(H_W)$, $H_W=\gamma_5 D_W$. $H_W$ is nearest-neighbor in time, but
$\mathrm{sign}(H_W)$ is a NON-LOCAL (dense-in-$t$) function of it. Hence:
- there is NO ultralocal temporal transfer matrix, and the exact Noether energy density is NOT a single
  time-link object;
- the naive $\tfrac12(\psi_{t+1}-\psi_{t-1})$ is the LEADING term; the true density carries tails that
  decay with the overlap locality range (exponential, set by the $H_W$ gap).

## Route A -- Noether "cut" construction (exact conserved $T_{00}$)
Promote the shift to a local one: displace fields/links only for $t\ge t_0$. Then $\delta S$ localizes across
the cut at $t_0$, and
$$
T_{00}(t_0) \;=\; \frac{\delta S}{\delta(\text{shift at } t_0)} .
$$
Nearest-neighbor part $\to$ standard hopping energy; the $\mathrm{sign}(H_W)$ part $\to$ a sum over pairs
$(t_a < t_0 \le t_b)$ weighted by the kernel $\langle t_a|\mathrm{sign}(H_W)|t_b\rangle$ (the tails).

## Route B -- Free / spectral (exact dispersion; pins the plateau)
Free limit: diagonalize in temporal frequency $\omega$ and spatial mode $p$; $D_{\rm ov}(\omega,p)$ is explicit
(2x2 for the 2-component field). The single-particle energies solve the propagator-pole condition
$$
\det D_{\rm ov}\big(\omega = i\,a_t E,\; p\big) = 0 .
$$
This gives the EXACT free energy per mode -- it directly settles whether the free $T_{00}$ state sits at
$3/R=0.567$ or is lattice-shifted (our single-op $0.569$ vs block-Hankel $0.533$ discrepancy). The density
$T_{00}$ then follows from the dispersion (group-velocity / $\partial_\omega D_{\rm ov}$ structure).

## Recommendation
Start with Route B: our validation IS the free limit, so it pins both the exact energy and the correct
operator normalization; then do Route A for the interacting conserved density.

## ALIGNED result (first-order / phase-space action, NM 2026-09-16)
First-order Euclidean action $S_E = \int(-i\,a + H\,dt)$: the TEMPORAL hops are the symplectic term $a$
(they set the canonical pairing / transfer matrix), and the SPATIAL hops ARE the Hamiltonian. So the energy
density is a SINGLE-TIMESLICE spatial Wilson-hop bilinear -- NO temporal derivative. From NM + `qed3int_v3-4`:
$$
H(t) = \sum_{\langle x,y\rangle\ \rm spatial\ nn}
   \eta^\dagger_x(t)\, D_{W;xy}\, \xi_y(t) \;+\; \xi^\dagger_x(t)\, D_{W;xy}^\dagger\, \eta_y(t),
\qquad D_{W;xy} = -\kappa_{xy} P_{xy} U_{xy},
$$
with the Wilson projector $P_{xy} = \tfrac12(1 - e^a_{xy}(x)\sigma_a)\,\Omega_{xy}$ (`qed3int_v3-4` Eq B.2;
analytic $\kappa,\Omega,e_a$ in `qed3_v2-6.pdf`). $\mathcal{O}_{T_{00}}(t)=H(t)$ (uniform sum) is already the
rotational scalar $\ell=0$; the $A_x Y_{00}$ weight is the same $\ell=0$ projection (constant at L1).

Propagator MISMATCH (`qed3int_v3-4`, meson bilinears Eq 5.4): only the cross-contractions are nonzero,
$$
\langle \xi\,\eta^\dagger\rangle = D_{\rm ov}^{-1} \equiv G, \qquad
\langle \eta\,\xi^\dagger\rangle = D_{\rm ov}^{-\dagger} \equiv \tilde G = 1 - D_{\rm ov}^{-1}\ \text{(massless GW)},
$$
and $\langle\xi\xi^\dagger\rangle=\langle\eta\eta^\dagger\rangle=0$. Write $W$ = the spatial Wilson kernel matrix.
Connected two-point (single loop; $\mathcal{O}=\eta^\dagger W\xi + \xi^\dagger W^\dagger\eta$, cross terms vanish):
$$
C_{T_{00}}(t,t_0) = -\,\mathrm{Tr}\big[W\,G(t,t_0)\,W\,G(t_0,t)\big]
                   \;-\; \mathrm{Tr}\big[W^\dagger\,\tilde G(t,t_0)\,W^\dagger\,\tilde G(t_0,t)\big].
$$
This is EXACTLY the code's $\sigma^2$ loop structure $S_4 + \tilde S_4$ (PS/FS), with the on-site density
vertex $S=\mathrm{diag}(w)$ REPLACED by the two-site spatial Wilson kernel $W$; distillation blocks:
$G\to$ `AblkS` (tau), $\tilde G\to$ the `taugw` block. So $T_{00}$ is a "hopping meson", reusing
`fs_gevp_point.make_config` with a matrix vertex $W$ instead of the scalar $w$.

## Open questions -- RESOLVED with NM (2026-09-16)
1. $H$ = the FIRST LINE of Eq (IV.1) of the free-limit paper `qed3_v2-6.pdf` -- the spatial Wilson operator,
   INCLUDING the $\tfrac12 r$-diagonal ("subtracts 1/2"). My earlier $\bar\psi\sigma_3(\psi_{t+1}-\psi_{t-1})$
   was the TEMPORAL/symplectic term (Eq IV.4, 2nd line) -- NOT the energy. So $W$ = spatial Wilson kernel:
   $$
   W_{xy} = \tfrac12\,\kappa_{xy}\,(-r\,\mathbb{1} + e^a_{xy}(x)\sigma_a)\,U_{xy}\,\Omega_{xy}\quad(x\ne y),
   \qquad W_{xx} = \tfrac12\, r\!\!\sum_{y\ \rm nn}\!\!\kappa_{xy}\,\mathbb{1},
   $$
   with the $M_5$ term DROPPED (overlap projection mass, not physical $H$). This is `dirac_simp.h:269-270`
   with $M_5=0$. $r=1$, $\Omega_{xy}=\cos\tfrac{\omega_{xy}}{2}+i\sigma_3\sin\tfrac{\omega_{xy}}{2}$ (Eq III.1).
2. Weight is ALREADY in $\kappa_{xy}=2A_{xy}/(\bar a_s\ell_{xy})$ (Eq IV.2). Use the link sum $\sum_{y_1,y_2}$
   as-is; NO extra $A_x Y_{00}$. $\mathcal{O}_{T_{00}}(t)=H(t)$ (the whole spatial operator), $\ell=0$ by
   rotational invariance.
3. Build $W$ in Python from `/mnt/barracuda22/qed3/qed3/geometry/data` (omega_n*, links_n*, alpha_n*,
   dualtriangleareas / link volumes, ell), reproducing `dirac_simp.h` set_kappa + the hop assembly. OK'd.
4. $\tilde G = \langle\eta\xi^\dagger\rangle$ block: PENDING -- NM says `taugw` is likely a FURNISHED
   propagator, not the plain dagger. Checking exact definition with "Fin: Two-meson" before coding this leg.

## FINDING (2026-09-16): drop the Wilson term -> naive e.sigma energy vertex
The full spatial Wilson $W=C+B$ splits into the naive part $C=\tfrac12\kappa\,(e^a_{xy}\sigma_a)\Omega$ (antihermitian)
and the Wilson $r$-term $B=\tfrac12\kappa(-r\mathbb{1})\Omega + \text{diag}$ (hermitian, spin-SCALAR). Free-limit tests:
- full $W$ ($r{=}1$): effmass slides to ~0.39 ($\sigma$-dominated via $B$).
- $W$ minus diagonal only: cleaner $\sigma$, 0.378 (worse) -- the scalar is the WHOLE $B$, not the diagonal.
- $W\to C$ ($r{=}0$, e.sigma hop only): plateau 0.563 at $dt$=18-20 -> $3/R=0.567$. Antihermitian to $7\times10^{-16}$.

Physics: the Wilson term is $O(a^2)$ with a ZERO at the physical pole (paper qed3_v2-6 Sec IV A) = doubler-removal
artifact; being a spin-scalar it overlaps $\sigma$ ($\Delta=2$). The genuine energy density is the NAIVE part $C$,
which is EOM-identical to the direct temporal form $\bar\psi\sigma_3 D_t\psi$ (both -> 0.567). DECISION (NM): use
$W\to C$, i.e. drop the Wilson term entirely (give up the projector, keep $e\!\cdot\!\sigma$ only). Coded as
`t00_stress_ham_claude.py` with `R=0` (default); `t00_wilson_kernel_claude.build_W(r=0)`.

## Corrected operator + contraction (final)
$$
\mathcal{O}_{T_{00}}(t) = \eta^\dagger(t)\,W\,\xi(t) + \xi^\dagger(t)\,W^\dagger\,\eta(t),\qquad
C_{T_{00}}(t,t_0) = -\mathrm{Tr}[W G(t,t_0) W G(t_0,t)] - \mathrm{Tr}[W^\dagger \tilde G(t,t_0) W^\dagger \tilde G(t_0,t)],
$$
$G=\langle\xi\eta^\dagger\rangle=D_{\rm ov}^{-1}$ (= `AblkS`), $\tilde G=\langle\eta\xi^\dagger\rangle$ (taugw block, def pending Q4).
Reuses `fs_gevp_point.make_config`; only the vertex changes: on-site $\mathrm{diag}(w)\to$ two-site matrix $W$.

## FINAL RESULTS + STATUS (2026-09-17)
Operator SETTLED: T_00 = the NAIVE e.sigma spatial Wilson hop, r=0 (Wilson -r*1 term = O(a^2) spin-scalar
doubler artifact that overlaps sigma; DROPPED). `O_H = eta^H W xi + h.c.`, ell=0 = full link sum (no Y00 for
the hop). Free limit -> 3/R = 0.567 VALIDATED.
- O_T (temporal `eta^H sigma3 D_t xi + h.c.`) is EXACTLY orthogonal to O_H (and to all on-site sigma_k) in the
  free limit -- discriminant = OFF-SITE-hop vs ON-SITE (not spin/hermiticity/derivative); overlap EOM
  (D_ov psi=0) != naive Dirac, so they sit in orthogonal sectors. INTERACTING O_T is S/N-limited (~0.2 by dt6,
  temporal derivative on gauge noise) -> UNUSABLE; cannot decide meson-vs-tensor.
- INTERACTING Nf2 gsq1.0 at0.2: `t00_stress_ham_interacting_claude.py` (reads <ENS>/ckpoint_lat.k spatial links
  via primal_links_n{L}_claude.dat from t00_dump_links_claude.cu). **L1 a_t*m=0.459(9) -> Delta_T00=2*m/m_axial
  =2.74** (m_axial=0.3346 = d4 axial tp ell1). **L2 (pole-GEVP)**: `t00_ham_pole_gevp_claude.py` = 12 five-fold
  vertex point-operators (mask_star) + global H, 13-op GEVP; **plain off[0] rebase to a SINGLE state
  (NKEEP=1, reb1@8 T0=3)** -> fit[9,14] a_t*m=0.578(41) -> **Delta_T00(L2)=3.16(23)** (m_axial=0.3664). Any
  Hankel offset ladder destabilizes the 13-op matrix (rank<=Nv=24); use plain off[0].
- CONCLUSION (measured, not over-claimed): Delta_T00 L1 2.74 -> L2 3.16(23) converges toward the protected
  Delta=3 at finer lattice = the L1 shortfall is O(a^2). Caveat: L2 Nv=24 TRUNCATED (of 84) + config set still
  growing -> not final; complete Nv=84 L2 peram (Fin:Two-meson) would nail it. Joint fermionic+gluonic GEVP
  (glue = qed3-a7) needed for the fully-conserved Delta=3. O_H^2 four-fermion deferred (complicates contraction).

## WHY T_00 DOES NOT COUPLE TO sigma^2 (or to any single-insertion source) -- ANTI-HERMITICITY (2026-09-18)
NM question: we now understand (with Fin: Two-meson) why $\sigma_{PS}^2$ does not couple to $\sigma_{PS}$ via
"GW-anti-hermiticity"; is $T_{00}$'s non-coupling the SAME argument? Answer: **same selection-rule STRUCTURE (odd
number of anti-hermitian factors $\Rightarrow$ correlator $=0$), but a DIFFERENT anti-hermitian object** -- $T_{00}$
uses an anti-hermitian VERTEX, the $\sigma^2$ case uses the anti-hermitian GW PROPAGATOR.

**$T_{00}$ side (verified).** The $r=0$ energy-density kernel is ANTI-HERMITIAN:
$$
W^\dagger = -W ,\qquad \Vert W + W^\dagger\Vert / \Vert W\Vert = 7\times 10^{-16}\ \text{(free AND with random unitary links)} .
$$
This is just the anti-hermiticity of the Dirac KINETIC operator (spatial $\gamma\!\cdot\!D$); dropping the Wilson
$+r\cdot\mathbf 1$ diagonal -- the ONLY hermitian piece -- is exactly what leaves it purely anti-hermitian. Unitary
gauge links preserve it, so it survives interactions. Hence the sink elemental obeys $\Phi_{W^\dagger}=V^\dagger
W^\dagger V = -\Phi_W$. The stress-tensor operator is $\mathcal O_{T_{00}} = \eta^\dagger W\xi + \xi^\dagger
W^\dagger\eta$, and any correlator $\langle \mathcal O_{T_{00}}(t)\, O(0)\rangle$ is LINEAR in the sink vertex:
$$
\langle \mathcal O_{T_{00}}\, O\rangle = \mathcal T[\Phi_W] + \mathcal T[\Phi_{W^\dagger}]
 = \mathcal T[\Phi_W] + \mathcal T[-\Phi_W] = 0 \quad\text{for ANY } O .
$$
Verified interacting (Nf2 g1.0 L1, per config): $\langle T_{00}\,\sigma\rangle\sim 10^{-20}$, $\langle T_{00}\,
\sigma^2\rangle\sim 10^{-18}$. By contrast $\langle T_{00}\,T_{00}\rangle$ is QUADRATIC in $W$, so the two halves
carry $(-1)^2=+1$ and ADD -- nonzero (the $a_t m\approx0.459$ signal). So the rule is: **linear/odd in $W$ vanishes,
quadratic/even survives.**

**$\sigma^2$ side (Fin, GW -- confirm exact statement with her).** Same reality argument, but the anti-hermitian
object is the CONTACT-SUBTRACTED GW propagator. The GW relation $\tau + \tau^\dagger = \mathbf 1$ gives
$$
G \equiv \tau - \tfrac12\mathbf 1 ,\qquad G^\dagger = -G \quad\text{(anti-hermitian)} .
$$
The single-meson $\to$ two-meson amplitude $\langle\sigma^2\,\sigma\rangle$ is ODD in the number of these
anti-hermitian $G$'s (the $\sigma$ vertex $M=\mathrm{diag}(w_{00})$ is HERMITIAN, $M^\dagger=+M$, so it does NOT
supply the sign -- unlike $W$), so it cancels; the diagonal $\langle\sigma^2\,\sigma^2\rangle$ / $\langle\sigma\,
\sigma\rangle$ is even and survives.

**Unified statement.** Write each correlator as $\mathrm{Re}\,\mathrm{Tr}$ of a loop of building blocks, each block
either anti-hermitian ($W$; or $G=\tau-1/2$) or hermitian ($M=\mathrm{diag}(w_{00})$; the sink/source scalar vertex).
Conjugating the trace and using the (anti-)hermiticity flips the overall sign by $(-1)^{n_{\rm anti}}$, so
$\mathrm{Re}\,\mathrm{Tr}=0$ whenever the total number of anti-hermitian factors is ODD. $T_{00}$ supplies its odd
factor as the VERTEX $W$ (kinetic anti-herm, not GW); $\sigma^2\!-\!\sigma$ supplies it as an odd count of GW
PROPAGATORS $\tau-1/2$. Same principle, different carrier -- so they are "the same argument" only at the level of
anti-hermitian-factor parity, NOT literally the same anti-hermiticity (W's is the Dirac kinetic one, independent of
the Ginsparg-Wilson relation).

### Two symmetries, distinct jobs (Fin's refinement) + the T_00 special case
The general statement is about $\mathrm{Re}\,\mathrm{Tr}$: (i) T-symmetry (real action, ensemble invariant under
$U\to U^*$) makes the T-odd $\mathrm{Im}$ part average to zero, which is what LICENSES taking $\mathrm{Re}$; (ii)
anti-hermiticity with an ODD count flips $\mathrm{Re}\,\mathrm{Tr}\to-\mathrm{Re}\,\mathrm{Tr}$, killing it. For
Fin's $\langle\sigma\,\sigma^2\rangle=\mathrm{Tr}[\Phi M\Phi M\Phi M]$ (3 anti-herm $M$'s, hermitian vertices $\Phi$),
$T^*=-T$ so $T$ is purely IMAGINARY -- nonzero per config -- and one needs BOTH: anti-herm kills $\mathrm{Re}$,
T-symmetry drops $\mathrm{Im}$.
**T_00 is stronger:** because $\mathcal O_{T_{00}}$ is LINEAR in the single anti-herm vertex and $\Phi_{W^\dagger}=
-\Phi_W$, the two h.c. halves cancel as FULL COMPLEX numbers, verified per config:
$\mathcal T[\Phi_W]+\mathcal T[\Phi_{W^\dagger}]\sim 3\times10^{-21}$ (both $\mathrm{Re}$ AND $\mathrm{Im}$), while a
single half $\mathcal T[\Phi_W]\sim 10^{-4}$. So $\langle T_{00}O\rangle=0$ needs ONLY vertex anti-hermiticity -- no
Re-projection, no T-symmetry.

### Corollary (NM): what DOES couple to T_00 needs a SECOND anti-hermitian kernel
$\langle T_{00}\,O\rangle$ carries exactly ONE anti-herm vertex $W$ (odd) whenever $O$ is built from HERMITIAN
vertices (scalar densities $\sigma_{PS}$, $\sigma^2$, identity-type), so it VANISHES. To get a nonzero correlator
$O$ must ALSO carry an anti-hermitian kernel, so the two anti-herm vertices sit in the SAME trace and the count is
EVEN: $(-1)^2=+1$, they ADD. This is precisely why $\langle T_{00}T_{00}\rangle\ne0$ (two $W$'s, quadratic). So T_00
couples only to other CURRENT-/KINETIC-like operators (anti-hermitian, "hop"/derivative kernels) -- e.g. itself, the
conserved currents -- and is orthogonal to every scalar-density (hermitian-vertex) operator. NECESSARY, not
sufficient: rotation/parity/flavor selection still apply on top (an anti-herm kernel of the wrong spin still gives 0).

## REVISION 2026-09-21/22 -- READ THIS BEFORE THE "FINAL RESULTS" SECTION ABOVE
The $\ell=0$ results above ("0.567 validated", "Delta_T00 2.74 -> 3.16") are RETRACTED as statements about $T$:
$\int T_{00}=H$ is a conserved charge and creates nothing; the spin-2 primary $T$ is the $\ell=2$ multiplet at $3/R$.
The old $\ell=0$ plateau was the $(2,2)$ meson {1,1,1,1}. Current status and evidence:
- $\ell=2$ projection of the SAME operator $O_H$ ($W_{ij}\to W_{ij}Y_{2m}(\hat n_{\rm mid})$, real $Y$, stays anti-hermitian):
  `t00_ell2_claude.py`, plan `t00_ell2_impl_plan_claude.md`.
- FREE $\Delta_T=3$ CONFIRMED with the EXACT free $D_{\rm ov}$ (`t00_free_exact_claude.py`, `t00_free_exact_impl_plan_claude.md`):
  $2m_{\ell2}/m_\sigma$ -> 2.963 ($a_t=0.2$), 2.990 ($a_t=0.1$), $a_t^2\to0$: 2.999. At L3/L4 the $\ell=0$ $O_H$ correlator is
  $t$-independent (conserved charge). Only $O_H$ (spatial hop) was used; $O_T$ (temporal $\sigma_3D_t$) was not needed.
- Connected-only $O_H$ = flavor-ADJOINT spin-2 ($3+\gamma_{\rm adj}$); the singlet needs the disconnected loops
  (`t00_ell2_singlet_claude.py`; L1 Nf2 g1.0: 1-3% of connected). "$\langle T_{00}\sigma^2\rangle=0$ exact" holds for the
  CONNECTED part only (the tadpole leg of the $W^\dagger$ half is $1-\tau$, so the loops add).
- Gauge energy density is $T^G_{00}=(E^2-B^2)/2g^2$ (NOT the action density $E^2+B^2$).
