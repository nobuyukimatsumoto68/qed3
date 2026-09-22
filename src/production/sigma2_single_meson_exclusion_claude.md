# Exclusion: $\sigma^2$ does not couple to a single (scalar) meson

> **CORRECTION (2026-09-17, supersedes the FS-triangle sections below).** The whole "$\sigma_{FS}^2$ couples
> at O(a)/O(1)" / "two inequivalent constructions" story below is an ARTIFACT of using the stored
> `peram/tau_gw` as the FS leg. `tau_gw` $=(1-D_{\rm ov}^\dagger)D_{\rm ov}^{-1}$ uses the FORWARD inverse and
> carries a bare $D_{\rm ov}^\dagger$; the correct FS leg (adjoint propagator $\langle\eta\xi^H\rangle=
> D_{\rm ov}^{-\dagger}$, dressed by $\tilde S=-(1-D_{\rm ov}^\dagger)$) **collapses to plain $\tau$** by GW
> ($\tilde S D_{\rm ov}^{-\dagger}=1-D_{\rm ov}^{-\dagger}=D_{\rm ov}^{-1}=\tau$). Full derivation:
> `fs_furnishing_derivation_claude.md`. VERIFIED: with the correct $\tau$ leg the FS triangle is $R\sim10^{-6}=0$
> (vs $-0.0046$ with $-\tau_{gw}$). **So $\sigma_{FS}^2=\sigma_{PS}^2$ (GW), both protected, NEITHER couples to a
> single meson.** There are NOT two constructions -- one correct leg ($\tau$), one wrong leg ($\tau_{gw}$).
> The sections below marked "O(a)/O(1)/two constructions" are retained only as the record of the wrong path.
> What still stands: $\langle\sigma_{00}\sigma^2\rangle=0$, and PS$==$FS.

> **MECHANISM CORRECTION (2026-09-18).** The exclusion mechanism below is stated as "$\sigma_3$-hermiticity,
> measure-independent." That is the **free-only** statement. The general, interaction-surviving mechanism is
> **GW anti-hermiticity of the normal-ordered propagator** $M=D_{\rm ov}^{-1}-\tfrac12$ ($M^\dagger=-M$ by GW):
> a closed loop of $n$ $M$'s with hermitian vertices obeys $T^*=(-1)^n T$, so the odd 3-loop triangle is purely
> imaginary $\Rightarrow \mathrm{Re}\,T=0$. This holds **config-by-config interacting** (2+1D 2-component has NO
> chirality, so there is no $\sigma_3$/$\gamma_5$-hermiticity of the $D_{\rm ov}^{-\dagger}=\sigma_3 D_{\rm ov}^{-1}\sigma_3$
> form). The **free** case additionally kills the imaginary part site-by-site (that is the free $\sigma_3$-herm
> identity below). Full derivation: `gw_antiherm_exclusion_mechanism_claude.md`.

## Statement

In the free limit, the four-fermion $0^{++}$ operator $\sigma^2_{00}$ has **zero overlap with any state
created by the local single scalar $\sigma_{00}=\bar\psi\psi$** ($\Gamma=\mathbb 1$, $\ell=0$):
$$
\boxed{\;\langle\,\sigma_{00}(t)\;\sigma^2_{00}(0)\,\rangle \;=\; 0\quad\text{exactly, at every }t\;}
$$
**General mechanism (holds interacting):** the normal-ordered propagator $M=D_{\rm ov}^{-1}-\tfrac12$ is
anti-hermitian by GW ($M^\dagger=-M$); a closed loop of $n$ $M$'s with hermitian vertices obeys $T^*=(-1)^n T$,
so the **3-propagator (odd)** triangle is purely imaginary and $\mathrm{Re}\,T=0$ — the physical correlator
vanishes, config-by-config (see `gw_antiherm_exclusion_mechanism_claude.md`).

**Free-only strengthening ($\sigma_3$-hermiticity):** in the *free* theory $D\sigma_3=-\sigma_3 D$
(Eq. C.16 of `qed3_v2-6.pdf`) gives $\sigma_3 G\sigma_3=-G$, so the $\sigma_3$-selection $\mathrm{Tr}[\Gamma_{\rm even}GGG]=0$
holds as an **algebraic (spinor-trace) identity, independent of the spatial measure** — site-by-site before the
area sums (it kills the *imaginary* part too). This is NOT available interacting (2+1D 2-component = no chirality);
there the $\ell=0$/GW-anti-herm route (above) is what remains, and it does.

## Direct verification (L1 AND L2, single free config)

`ground_coupling_check_claude.py` measures the triangle $C_{12}=\langle\sigma_{00}(t)\,\sigma^2_{00}(0)\rangle$
and the cross $\langle\sigma_{00}\,O_A\rangle$, normalized by $\langle\sigma_{00}\sigma_{00}\rangle$:

| | $C_{12}$ (triangle) | $C_{12}/\langle\sigma_{00}\sigma_{00}\rangle$ | $\langle\sigma_{00}O_A\rangle/\langle\sigma_{00}\sigma_{00}\rangle$ |
|---|---|---|---|
| **L1** (Nv24, uniform areas) | $\sim10^{-10}$-$10^{-12}$ (noise) | $\sim10^{-6}$ | $\sim10^{-5}$ |
| **L2** (Nv84, non-uniform areas) | $\sim10^{-12}$-$10^{-13}$ (noise) | $\sim10^{-6}$ | $\sim10^{-5}$ |

Both vanish to numerical noise, at **both** L, and the effmass of $C_{12}$ is pure noise (no plateau). So
$\langle0|\sigma^2_{00}|{\rm single~scalar~meson}\rangle=0$ at both finite lattices. **The exclusion is exact
and measure-independent** -- it does NOT break at L2.

## What this corrects (the "ground leak" misreading)

Earlier the per-diagram diagram-A effmass (`per_diagram_ABE_free_claude.py`) was seen to plateau at
$2E_1=0.556$ at L1 (the $(2,2)$) but to **slide** at L2 through $(2,2)$ down to $\approx0.39$. That $0.39\approx
m_\sigma=2E_0$ was mis-read (by me, and provisionally in the codeset audit) as diagram A coupling to the single
scalar **ground** meson, with the L1 cleanliness dismissed as a "uniform-measure accident."

**That interpretation is wrong.** The direct triangle test above shows $\sigma^2$ has NO overlap with the
single scalar meson at either L. Since $\sigma_{00}$ does not overlap the $0.39$ object (triangle $=0$), the
$0.39$ in the diagram-A DIAGONAL is **not** the scalar ground. Two candidates remain, both consistent with the
vanishing triangle:
1. **Noise floor** -- the L2 diagram-A correlator flattens to a constant $\sim2\times10^{-13}$
   (`diagramA_corr_log_L2`); the effmass sliding through $0.39$ is heading to $0$, not a plateau.
2. **$\ell$-mixing artifact** -- the non-uniform L2 measure leaks $\ell=1$ into the "$\ell=0$" $\sigma^2_{00}$;
   the $\sigma_3$-odd diagram A can then pick up the $\ell=1$, $\sigma_3$-odd state that is energy-degenerate
   with $2E_0=0.393$ but is a DIFFERENT state (correctly not overlapped by the $\ell=0$ scalar $\sigma_{00}$).

Either way it is not a "two-meson couples to a single scalar meson." That coupling is zero at every lattice.

## Scope: what IS coupled

The exclusion is specifically about the $\sigma_3$-even scalar tower that $\sigma_{00}$ creates. The physical
single-meson content of $\sigma^2$ -- the $(2,2)=\{1,1,1,1\}$ excitation ($\Delta=4$), reached through
diagram A's $\sigma_3$-**odd** nonlocal kernel $\tilde\tau$ (the $O_A$ construction) -- is a separate channel
and is NOT excluded; $\langle O_A O_A\rangle\to2E_1=0.556$ at L1. (Why $O_A$ shows $2E_1$ and no $2E_0$ despite
$\tilde\tau$'s even $-\tfrac12$ part remains the open puzzle of `sigma_quantum_numbers_claude.md:126`; the
triangle result here says only that whatever $O_A$/diagram A couples to, it is orthogonal to the local
$\sigma_{00}$ scalar tower.)

## $\text{PS}^2 \to \text{FS}$ (measured) -- nonzero, but an $O(a)$ parity-impurity artifact

The exclusion above uses a $\sigma_3$-**even** sink vertex ($\sigma_{PS}$). The FS furnishing $\tilde S=-(1-
D_{\rm ov}^\dagger)$ carries a $\sigma_3$-**odd** part, so the triangle $\langle\sigma_{FS}(t)\,\sigma_{PS}^2(0)
\rangle$ is NOT protected by $\mathrm{Tr}[\Gamma_{\rm even}GGG]=0$. Measured with Fin's furnishing recipe
(`fs_ps2_triangle_recipe_claude.md`; only the sink-arriving leg is furnished $\to-\tau_{gw}$, prefactor $-1$):

| | validation $\tau_{gw}\!\to\!\tau$ (must $\to0$) | $R=\langle\sigma_{FS}\sigma_{PS}^2\rangle/\langle\sigma_{00}\sigma_{00}\rangle$ ($dt\!\sim\!10$) |
|---|---|---|
| L1 (uniform) | $\sim10^{-11}$ ✓ | $\approx-0.0046$ |
| L2 (non-uniform) | $\sim10^{-12}$ ✓ | $\approx-0.00026$ |

So $\sigma_{PS}^2$ **does** couple to a single $\sigma_{FS}$ (nonzero, unlike PS$\to$PS) -- the FS furnishing's
$\sigma_3$-odd part breaks the protection. But the coupling is $\sim18\times$ **smaller at L2 than L1**: it
shrinks strongly as the lattice fines, exactly the signature of an **$O(a)$ parity impurity** of $\sigma_{FS}$
($\sigma_{FS}$ has no definite lattice parity; it sharpens only as $a\to0$), extrapolating to $0$ in the
continuum (parity restoration). Driver `ps2_fs_triangle_claude.py` (validation check built in).

## Full flavor table: all four $\langle\sigma_X(t)\,\sigma_Y^2(0)\rangle$ triangles

Measured L1 and L2 (`flavor_triangle_duality_claude.py`), ABSOLUTE numerator effmass + amplitude
$R=C/\langle\sigma_{00}\sigma_{00}\rangle$ (the relative single-meson overlap; amplitude-vs-$L$ is the
$O(a)$-vs-$O(1)$ decider). Furnishing rule: a leg is $-\tau_{gw}$ iff its arrival (row) vertex is FS; the FS
equal-time **source** leg is $-\tau_{gw}(s,s)+c_{FS}I$ with $c_{FS}=\mathrm{tr}[\Phi\tau_{gw}]/\mathrm{tr}[\Phi]$
($=-0.2256$ L1) -- there is NO $\pm\tfrac12$ analogue for FS (sign matters: use $-\tau_{gw}$ CONSISTENTLY;
the tadpole test $\mathrm{tr}[\Phi\,\tilde\tau^{FS}]=0$ does NOT catch the sign since $0=-0$).

| triangle (sink$\leftarrow$source) | effmass | $R$(L1) | $R$(L2) | verdict |
|---|---|---|---|---|
| $\sigma_{PS}\leftarrow\sigma_{PS}^2$ | noise | $\sim10^{-6}$ | $\sim10^{-6}$ | **exactly 0** ($\sigma_3$-herm) |
| $\sigma_{FS}\leftarrow\sigma_{PS}^2$ | $\to m_\sigma$ region | $-0.0031$ | $-0.00005$ | **$O(a)$** ($\sim60\times$ shrink) |
| $\sigma_{PS}\leftarrow\sigma_{FS}^2$ | $\to m_\sigma$ | $+0.052$ | $+0.034$ | **$O(1)$** ($\sim1.5\times$) |
| $\sigma_{FS}\leftarrow\sigma_{FS}^2$ | $\to m_\sigma$ | $-0.037$ | $-0.033$ | **$O(1)$** ($\sim1.1\times$) |

Physical constraint satisfied everywhere: no channel is lighter than $m_\sigma$ and none reaches the two-meson
$0.756$ (a 2-fermion sink cannot overlap a 4-fermion state) -- so the plateaus are real single-meson states.
**No "FS$\leftarrow$FS $=0$" duality** (the PP=FF even-loop equality is a four-point property, does not extend
to these single-meson triangles).

### DECISIVE test -- is the FS-source $O(1)$ genuine or the unsolved equal-time contact?

Time-split the FS source (src1 @ $0$, src2 @ $\Delta$): the src-src leg becomes a time-separated
$-\tau_{gw}(0,\Delta)$ with NO equal-time contact (`ps_from_fs_timesplit_claude.py`). Result (L1, PS$\leftarrow$FS):
$\Delta=1$ $R\!\sim\!0.086$, $\Delta=2$ $R\!\sim\!0.032$, both effmass $\to m_\sigma$ -- the coupling
**PERSISTS** (does not die for $\Delta>0$). So the $O(1)$ is the $\sigma_3$-odd **furnishing structure**
(present on every leg), NOT the equal-time contact. **Genuine.**

## THE HEADLINE: two INEQUIVALENT $\sigma_{FS}^2$ constructions (Fin, reconciling the PP==FF tension)

The $O(1)$ result above is in direct tension with the early **PP==FF bit-identical** finding
(`sigma2_flavorgeom_FULL_free` cache: $C[0,0]=C[3,3]$ to $\sim10^{-18}$). If explicit-$\tau_{gw}$ $\sigma_{FS}^2$
couples to the single meson at $O(1)$, then $FF=\langle\sigma_{FS}^2\sigma_{FS}^2\rangle$ has a single-meson pole
$|\langle{\rm meson}|\sigma_{FS}^2\rangle|^2\neq0$ -- but $PP$ (protected) does not, so $PP$ and $FF$ could NOT
be bit-identical. The resolution: **they are different objects.**

- **flavfac-$FF$ (the `sigma2_flavorgeom` CACHE):** built AblkS-only $+$ per-loop factor
  $(1+(-1)^{n_{FS}})$ -- **plain $\tau$ legs, NO $\tau_{gw}$.** This has NO $\sigma_3$-odd furnishing content,
  so it does NOT couple to the single meson -> flavfac-$FF ==\,PP$ (both protected). *That* is why PP==FF is
  bit-identical.
- **explicit-$\tau_{gw}$ $\sigma_{FS}^2$ (the Eq 5.1-5.3 definition, = production `four_point`'s
  $S_4[-\tau_{gw}]$):** carries the $O(1)$ $\sigma_3$-odd content and DOES couple to the single meson.

So **flavfac-$FF$ $\neq$ explicit-$\tau_{gw}$ $\sigma_{FS}^2$** -- the flavfac collapse silently DROPS the $O(1)$
$\sigma_3$-odd furnishing content. This also explains the very first discrepancy of the whole investigation:
flavfac $\sigma^2_{00}$ $\to0.56$ (protected, $(2,2)$) vs explicit `fs_gevp_point` $\to0.378$ (couples to the
single meson). "$\sigma_{FS}^2$" is **ambiguous in the codebase: two inequivalent constructions.**

**IMPLICATION (flag to NM + the $\sigma^2$-$F^2$ mixing thread):** the flavorgeom-cache $FF$ channel (and the
6x6 flavor$\times$geometry GEVP built on it) is the **flavfac/protected** version, NOT the Eq-5.1-5.3 FS
operator. NM should decide which the cache SHOULD implement. (Co-check available: does production `four_point`
FS.FS effmass show $m_\sigma$? -> tells which construction production actually uses where.)

## Complete picture (FINAL)

- $\sigma_{PS}^2 \to$ single $\sigma_{PS}$: **exactly zero at all $a$** ($\sigma_3$-hermiticity, protected).
  $\sigma_{PS}^2$ is the clean $(\bar\psi\psi)^2$ interpolator.
- $\sigma_{PS}^2 \to$ single $\sigma_{FS}$: nonzero but **$O(a)$** parity impurity ($\to0$ in continuum).
- **explicit-$\tau_{gw}$** $\sigma_{FS}^2 \to$ single $\sigma_{PS}$ and $\to\sigma_{FS}$: **genuine $O(1)$**
  coupling to the single scalar meson at $m_\sigma$, via the $\sigma_3$-odd FS furnishing $(1-D_{ov}^\dagger)$
  (time-split confirmed: furnishing, not equal-time contact). Does NOT vanish in the continuum.
- **flavfac-$FF$** (cache): $==PP$, protected, NO single-meson coupling. (Inequivalent to the above.)

**Consequence (Chester-Pufu $(\bar\psi\psi)^2$ / GEVP basis):** $\sigma_{PS}^2$ is the clean $(\bar\psi\psi)^2$
(protected either way). The explicit-$\tau_{gw}$ $\sigma_{FS}^2 = (\bar\psi\psi)^2 +$ an $O(1)$ furnishing
(Dirac-op/derivative-like) piece that couples to the single meson -- do NOT use it in a basis meant to isolate
$(\bar\psi\psi)^2$. The flavfac-$FF$ happens to also be protected, but it is NOT the Eq-5.1-5.3 FS operator.

## Aside: $\langle\sigma_{PS}(t)\,\sigma_{FS}(0)\rangle$ single-meson cross (NM asked; thought it was 0)

NOT zero. `sigma_ps_fs_2pt_claude.py` (L1): $R_{PF}=C_{PF}/C_{PP}\to\sim-0.8$ (toward $-1$) -- $\sigma_{PS}$ and
$\sigma_{FS}$ overlap the SAME single meson with opposite sign. Reason: $\Gamma_{FS}$ = $\sigma_3$-even
($\sigma_{00}$) $+$ $\sigma_3$-odd (furnishing); the even$\times$even piece is the ordinary $\sigma$ two-point
(nonzero, $\to m_\sigma$), the even$\times$odd piece vanishes by $\mathrm{Tr}[\Gamma_{\rm even}G\Gamma_{\rm odd}G]=0$.
So the odd furnishing drops from the CROSS and the ordinary meson overlap survives. (What IS zero: the
$\sigma_3$-odd part alone.)

## Files / provenance

- `ground_coupling_check_claude.py` -- triangle $C_{12}$ + cross $\langle\sigma_{00}O_A\rangle$, L1 & L2.
- `per_diagram_ABE_free_claude.py`, `diagramA_corr_linear_claude.py` -- the diagram-A diagonal that motivated
  (and whose $0.39$ was mis-read).
- Selection rule: $\sigma_3$-hermiticity, `o_a_operator_note_claude.md:37-42`; open puzzle
  `sigma_quantum_numbers_claude.md:126`.
