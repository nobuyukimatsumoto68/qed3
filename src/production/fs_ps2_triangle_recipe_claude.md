# $\langle \sigma_{FS}(t)\,\sigma_{PS}^2(0)\rangle$ triangle -- furnishing recipe

> **SUPERSEDED / CORRECTED (2026-09-17).** Two errors in this note, both fixed in
> `fs_gw_collapse_v_agent_two_meson_claude.md` (read that instead):
> 1. **No $\sigma_3$-hermiticity.** This system does NOT have $\sigma_3$- (gamma-) hermiticity (3D / 2-component,
>    no chirality). Every statement below that attributes a zero to "$\sigma_3$-hermiticity" or "$\mathrm{Tr}[\Gamma_\text{even}GGG]=0$"
>    or calls a leg "$\sigma_3$-odd/even" is WRONG. The suppression $\langle\sigma^2|m_{PS}\rangle\approx0$ is empirical
>    at the COMPLETE distillation basis (and leaks under truncation), mechanism not derived.
> 2. **The $O(a)/O(1)$ FS-coupling conclusions are RETRACTED.** The apparent $\sigma_{FS}^2\to$ single-meson coupling
>    was a `tau_gw` artifact (the FS leg was built from the FORWARD inverse, carrying a spurious bare $D_{ov}^\dagger$).
>    The correct FS leg collapses to plain $\tau$ by GW, so $\sigma_{FS}==\sigma_{PS}$ and there is NO genuine FS
>    single-meson coupling. See the companion writeup.
>
> What survives: the S/S~ split (Eq 5.5) and the observation that `tau_gw` (forward-furnished) is the wrong leg.

Requested by the "{1,1,1,1}" agent (NM: "PS$^2$ to FS is worth checking"). Single $\sigma_{FS}$ sink
bilinear, two-$\sigma_{PS}$ source. Goal: the block expression with the correct FS furnishing, analogous to
their PS triangle
$$
C_{12}=\langle \sigma_{00}(t)\,\sigma^2_{00}(0)\rangle = -2\,\mathrm{Tr}[\Phi_0\,\tilde\tau_{00}\,\Phi_0\,\tau_{st}\,\Phi_t\,\tau_{ts}] = 0 .
$$

## 1. The S / S~ split (qed3int_v3-4.pdf Eq 5.5)

$\sigma = \eta^\dagger S\xi + \xi^\dagger \tilde S\eta$, with $S=1$, $\tilde S=1$ (PS) / $-(1-D_{ov}^\dagger)$ (FS).
A single fermion loop through $k$ bilinears must use the SAME vertex type ($S$ or $\tilde S$) at every vertex --
$S$ and $\tilde S$ never mix on one loop. So the three-point splits exactly like the four-point:
$$
\langle \sigma_{FS}\,\sigma_{PS}\,\sigma_{PS}\rangle = \langle S_3\rangle + \langle \tilde S_3\rangle .
$$

- $\langle S_3\rangle$: every vertex uses $S=1$ (the FS sink uses its plain $\eta^\dagger\xi$ term). ALL legs plain
  $\tau$. This is exactly the plain triangle $=\tfrac12 C_{12}=0$ (protected by $\sigma_3$-hermiticity,
  $\mathrm{Tr}[\Gamma_\text{even}GGG]=0$, exact and measure-independent).
- $\langle \tilde S_3\rangle$: every vertex uses $\tilde S$. Sources are PS so $\tilde S=1$ (plain $\tau$); the
  sink is FS so $\tilde S=-(1-D_{ov}^\dagger)$ (furnished). This is NOT protected -- $(1-D_{ov}^\dagger)$ is
  $\sigma_3$-MIXED ($\sigma_3(1-D_{ov}^\dagger)\sigma_3 = 1-D_{ov}$), so the even-vertex zero does not apply.

Therefore
$$
\langle \sigma_{FS}(t)\,\sigma_{PS}^2(0)\rangle = \langle \tilde S_3\rangle \neq 0 \quad(\text{generically}) ,
$$
confirming the agent's expectation.

## 2. Furnishing rule (matches the production FS convention)

The furnished perambulator (chunk 2) is
$$
\tau_{gw}(a,b) = V^\dagger(a)\,(1-D_{ov}^\dagger)\,D_{ov}^{-1}\,V(b) ,
$$
i.e. $(1-D_{ov}^\dagger)$ is applied at the ROW (sink) time $a$. In `compute_diags` the leg object is
`legs[a_snk, a_src]`, and the production four-point furnishes FS with `legs = -taugw` on ALL legs
(`distill_contract_claude.py:284`, `four_point()` FS.FS $= S_4[\tau] + S_4[-\tau_{gw}]$).

RULE (equivalent to "all legs $-\tau_{gw}$ when all vertices FS"): a leg `legs[a_snk, a_src]` is furnished
($\to -\tau_{gw}(a_\text{snk}, a_\text{src})$) iff the vertex at its ARRIVAL (row) endpoint $a_\text{snk}$ is FS;
otherwise plain $\tau$.

For the mixed triangle: only the leg ARRIVING at the FS sink (row $=t$) is furnished; every leg arriving at a
PS source (row $=0$) stays plain.

## 3. The block expression

Loop: src1 $\xrightarrow{\tilde\tau_{00}}$ src2 $\xrightarrow{\text{(arrives sink)}}$ sink
$\xrightarrow{\text{(arrives src1)}}$ src1. Only the src2$\to$sink leg arrives at the FS sink.
$$
\boxed{\;\langle \sigma_{FS}(t)\,\sigma_{PS}^2(0)\rangle
= -\,\mathrm{Tr}\!\big[\Phi_0\,\tilde\tau_{00}\,\Phi_0\,\tau_{st}\,\Phi_t\,(-\tau_{gw})_{ts}\big]\;}
$$
with (in `compute_diags` leg names, source $a_\text{src}=0$, sink $a_\text{snk}=dt$):

- $\Phi_0,\Phi_t$ = the $Y_{00}$ area vertex $V^\dagger \mathrm{diag}(A_x Y_{00}) V$ at source / sink.
- $\tilde\tau_{00} = \tau(0,0) - \tfrac12 I$ -- the PS source--source equal-time leg (plain GW contact $-\tfrac12$).
  Note the sink is a SINGLE bilinear: it has NO equal-time self-leg, so the FS equal-time contact
  $-\tfrac12(\tau_{gw}+\tau)$ does NOT appear here.
- $\tau_{st} = $ `legs[0, dt]` (arrives at the source) -- plain $\tau$.
- $(-\tau_{gw})_{ts} = -\tau_{gw}(dt, 0) = $ `-taugw[dt, 0]` (arrives at the FS sink) -- the ONE furnished leg.

Prefactor $-1$ (single loop), NOT $-2$: the PS $-2$ was $\langle S_3\rangle + \langle\tilde S_3\rangle$ with the
two equal; here $\langle S_3\rangle=0$ and only $\langle\tilde S_3\rangle$ (one term) survives.

## 4. Checks before trusting a number

1. Set `taugw -> tau` in the boxed expression: it must reproduce $\tfrac12 C_{12}=0$ (all-plain, protected). A
   nonzero result there = a leg/sign bug.
2. Validate the $\tau_{gw}$ SIGN/orientation against production `four_point()`: your all-FS analogue (furnish
   BOTH the sink-arriving AND source-arriving legs) must match the FS.FS building blocks that already
   reproduce the validated $C_S^{FS}$ effmass $0.32$.
3. The result should be REAL up to the usual FS non-Hermiticity; take Re.

## 5. Physics expectation

$\langle\sigma_{FS}\sigma_{PS}^2\rangle$ is nonzero only through the $\sigma_3$-ODD part of $(1-D_{ov}^\dagger)$,
which is the lattice parity impurity of $\sigma_{FS}$ ("$\sigma_{FS}$ has no definite parity on the lattice;
sharpens only as $a\to0$", `project_scalar_ylm_corr`). So this triangle is an $O(a)$ PARITY-IMPURITY probe: it
should be nonzero at finite $a$ and extrapolate to $0$ as $a\to0$ (parity restoration), and it should shrink as
the mesh regularizes (L increases / $a_t$ decreases). If instead it is $a$-independent, that flags a genuine
parity-odd content, not an artifact -- worth reporting either way.

Refs: qed3int_v3-4.pdf Eq 5.1-5.5; `distill_contract_claude.py` `compute_diags`/`four_point`;
`project_scalar_ylm_corr` (FS parity mixing); the PS triangle exclusion `sigma2_single_meson_exclusion_claude.md`.

## 6. FS SOURCE ($\sigma_{FS}^2$) and the equal-time contact (the growing-R bug)

The FS SOURCE has an equal-time source--source leg (the two bilinears at $t=0$). Unlike the PS case
($\tilde\tau=\tau-\tfrac12$, tadpole exactly zero), the FS equal-time leg carries a NONZERO tadpole. Measured
on the free perams ($\Phi=V^\dagger\mathrm{diag}(A Y_{00})V$, $\mathrm{tr}\,\Phi=7.0898$):

| quantity | L1 (Nv24) | L2 (Nv84) |
|---|---|---|
| $\mathrm{tr}[\Phi\,\tau_{00}]$ (PS raw) | $3.5449=\tfrac12\mathrm{tr}\Phi$ | $3.5449$ |
| $\mathrm{tr}[\Phi(\tau_{00}-\tfrac12)]$ (PS sub) | $\sim10^{-10}$ (=0) | $\sim10^{-9}$ (=0) |
| $\mathrm{tr}[\Phi\,\tau_{gw,00}]$ (FS raw) | $-1.5994$ | $-1.7582$ |
| $\mathrm{tr}[\Phi(-\tau_{gw,00})]$ (leg $=-\tau_{gw}$) | $+1.5994$ (NONzero) | $+1.7582$ |
| $\mathrm{tr}[\Phi(-\tfrac12(\tau_{gw}+\tau))]$ (AblkSt) | $-0.9727$ (NONzero) | $-0.8933$ |
| $c_{FS}=\mathrm{tr}[\Phi\tau_{gw,00}]/\mathrm{tr}\Phi$ | $-0.22560$ | $-0.24799$ |

- The PS $-\tfrac12$ works because $\mathrm{tr}[\Phi\tau_{00}]=\tfrac12\mathrm{tr}\Phi$ EXACTLY (GW). There is NO
  $\pm\tfrac12$ analogue for FS: both candidate FS equal-time legs ($-\tau_{gw}$ and the fs_gevp_point
  $-\tfrac12(\tau_{gw}+\tau)$) leave an $O(1)$ tadpole. That un-subtracted source tadpole is the spurious
  slowly-varying / constant piece that makes $R$ GROW with $dt$ in PS$\leftarrow$FS and FS$\leftarrow$FS.
- The contact that zeroes the $\Phi$-weighted FS tadpole is $c_{FS}=\mathrm{tr}[\Phi\tau_{gw,00}]/\mathrm{tr}\Phi
  = -0.2256$ (L1), $-0.2480$ (L2) -- L-DEPENDENT, not universal. This IS the documented "FS forward-$\tau'$
  contact $\neq\tfrac12$, UNSOLVED" issue, now quantified.

FIX: vacuum/one-point-subtract the FS source. Either (i) use the FS equal-time source leg
$\tau_{gw}(0,0) - c_{FS}\,I$ with $c_{FS}=\mathrm{tr}[\Phi\tau_{gw}(0,0)]/\mathrm{tr}\Phi$ measured PER CONFIG,
or (ii) subtract the connected triangle's $dt\to\infty$ constant (correlator-level vacuum subtraction). The
source-arriving-leg furnishing ($-\tau_{gw}$, arrival vertex FS) was correct; ONLY the equal-time contact was
wrong. After the fix, $R$ should decay (no growing constant).

## 7. Construction choice (Q1) and the duality (Q2)

Q1 -- use EXPLICIT $\tau_{gw}$ furnishing, NOT flavfac. The flavfac form $\prod_\text{loops}(1+(-1)^{n_{FS}})$
is a $\sigma^2\times\sigma^2$ (closed 4-point) collapse. Applied to a single-meson triangle (one 3-vertex loop)
it gives $(1+(-1)^{n_{FS}})$ with $n_{FS}=1$ for FS$\leftarrow$PS $\Rightarrow 0$ -- which CONTRADICTS the
certified FS$\leftarrow$PS$=-0.0046$. So flavfac does not govern these triangles; keep explicit $\tau_{gw}$.

Q2 -- the PP$=$FF / even-loop PS$=$FS duality NM recalls is a property of the CLOSED $\sigma^2\times\sigma^2$
four-point (even cycles), NOT the single-meson triangles. It does NOT predict FS$\leftarrow$FS$=0$. The
triangles are governed by explicit furnishing + parity:
- PS$\leftarrow$PS $=0$ exactly ($\sigma_3$-herm, measure-independent).
- FS$\leftarrow$PS, PS$\leftarrow$FS, FS$\leftarrow$FS: $O(a)$ parity impurities (the $\sigma_3$-odd part of
  $(1-D_{ov}^\dagger)$), each $\to0$ as $a\to0$. Whether FS$\leftarrow$FS is extra-suppressed ($O(a^2)$, two FS
  insertions) or same order is EMPIRICAL -- test after the source fix. Do not assume FS$\leftarrow$FS$=0$ from
  the four-point duality.
