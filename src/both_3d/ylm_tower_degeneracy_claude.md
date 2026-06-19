# Per-m degeneracy of the local-current $Y_{\ell m}$ tower: rep-theory, not a bug

## TL;DR
The $Y_{\ell m}$ tower code is **correct**. The per-m non-degeneracy we saw is genuine representation
theory, in two parts:
1. **Scalar channel** $s_3 = J\!\cdot\!\hat n$ (normal/temporal current) is a true scalar field, so its
   scalar-$Y_{\ell m}$ tower degenerates exactly according to the icosahedral irrep content of $D^\ell$:
   $\ell{=}1\to T_1$ (degenerate), $\ell{=}2\to H$ (degenerate), $\ell{=}3\to T_2\oplus G$ (**reducible -> splits**).
2. **Tangent channels** $s_1{=}\hat\theta,\ s_2{=}\hat\phi$ are components of a tangent **vector field**.
   The code projects them onto **scalar** $Y_{\ell m}$, which is not their irrep basis, so $s_1,s_2$ and
   $sp=s_1s_1+s_2s_2$ are not icosahedral irreps and are **not** per-m degenerate -- even at $\ell=1$, even
   at zero noise. The covariant decomposition of a tangent vector field uses **vector spherical harmonics**.

No code change is required for correctness; only a *re-interpretation*. If a per-m-degenerate spatial tower
is wanted, project the tangent current onto vector spherical harmonics (Sec. 5).

## 1. What was verified correct
- **`Ylm_real`** (includes/valence_claude.h): passes the addition theorem
  $\sum_m Y_{\ell m}(\hat a)Y_{\ell m}(\hat b)=\frac{2\ell+1}{4\pi}P_\ell(\hat a\!\cdot\!\hat b)$ and the
  self-sum $\sum_m Y_{\ell m}^2=\frac{2\ell+1}{4\pi}$ to machine precision ($\le 3\times10^{-17}$) for
  $\ell=0\!-\!3$. So the real harmonics, per-m normalization, and `mult_Ylm_real` area weighting are right.
- The estimator, area weights $A_n=$`dual_areas`, and per-m bookkeeping are confirmed by Sec. 2 below.

## 2. Scalar channel $s_3=J\!\cdot\!\hat n$ reproduces icosahedral irreps exactly
On $S^2\subset\mathbb{R}^3$ the **normal** direction at site $n$ is $\hat n$ itself (globally, covariantly
defined), so $J\!\cdot\!\hat n$ is a genuine **scalar field**. (Equivalently $s_3$ is the temporal current
on the cylinder, with the global time axis.) Its scalar-$Y_{\ell m}$ projection
$\Sigma^3_{\ell m}=\sum_n A_n Y_{\ell m}(\hat n)\,(J\!\cdot\!\hat n)(n)$ therefore transforms as a pure $D^\ell$,
and the per-m two-point must degenerate according to the icosahedral content of $D^\ell$:

$$D^0=A,\quad D^1=T_1,\quad D^2=H,\quad D^3=T_2\oplus G,\quad D^4=G\oplus H,\ \dots$$

Deterministic free (U=1, exactly icosahedral), $s_3$ per-m at $dt=5$:

| $\ell$ | icos | per-m values ($\times10^{-5}$) | degeneracy |
| --- | --- | --- | --- |
| 1 | $T_1$ | $-4.378,\,-4.378,\,-4.378$ | all equal (3) |
| 2 | $H$ | five identical $=-3.461$ | all equal (5) |
| 3 | $T_2\oplus G$ | $-2.84,\,-4.26,\,{\sim}0,\,-7.11,\,{\sim}0,\,-4.26,\,-2.84$ | **splits** (not 7-fold) |

$\ell{=}1,2$ are fully degenerate; $\ell{=}3$ splits, exactly because $D^3$ is **reducible** under the
icosahedral group. This is a positive confirmation that the harmonic projection + areas + per-m machinery
are correct (a bug would not reproduce the irrep pattern).

**Takeaway 1 (a likely "misunderstanding"):** full per-m degeneracy is only guaranteed while $D^\ell$ is
icos-irreducible, i.e. $\ell\le 2$. From $\ell=3$ on, even the scalar tower splits into sub-multiplets.

## 3. Tangent channels are a vector field, not a scalar
The spatial current is a **tangent vector field** $J_{\rm tan}=J_\theta\,\hat\theta+J_\phi\,\hat\phi$.
The code (correctly, in the orthonormal vierbein) stores $s_1=J_\theta$, $s_2=J_\phi$, then projects each
*frame component* onto the **scalar** harmonic:
$$\Sigma^{1}_{\ell m}=\sum_n A_n Y_{\ell m}(\hat n)\,J_\theta(n),\qquad
  \Sigma^{2}_{\ell m}=\sum_n A_n Y_{\ell m}(\hat n)\,J_\phi(n).$$
A scalar harmonic is the right basis only for a **scalar** field. For a vector field the rotation group
acts on the frame index too: under an isometry $R$, $J^a(n)\to \Lambda(R,n)_{ab}\,J^b(Rn)$ with a
**position-dependent** frame rotation $\Lambda$. Hence $\Sigma^{a}_{\ell m}$ does not transform as a pure
$D^\ell$, and the per-m two-point is not an icos irrep.

Two consequences, both seen in the exact free data:
- **$s_1,s_2$ (and $sp$) are not per-m degenerate even at $\ell=1$.** Deterministic free, $sp$ per-m,
  rel-spread grows with $dt$ (12% at $dt{=}3$ -> 50% at $dt{=}8$), with *different per-m effective masses*
  $m_{\rm eff}(dt{\approx}8)\approx 0.47,0.50,0.54$ for $m={-}1,0,{+}1$ -- different IR states, not amplitudes.
- **$g^{s_1}_\ell \ne g^{s_2}_\ell$** (m-summed): $\hat\theta\hat\theta\ne\hat\phi\hat\phi$, e.g. $\ell{=}1$,
  $dt{=}10$: $|g^{s_1}|{=}1.7,\ |g^{s_2}|{=}3.9,\ |g^{s_3}|{=}3.9\,(\times10^{-6})$. This $\hat\theta/\hat\phi$
  anisotropy is the coordinate-frame artifact of projecting a vector field on a scalar harmonic.

**Takeaway 2 (the main "misunderstanding"):** $sp=s_1s_1+s_2s_2$ built from scalar-$Y_{\ell m}$-projected
frame components is **not** a rotational scalar (the same-frame-index contraction between two different
points is not covariant -- the frames there are not parallel transported). So per-m $sp$ should *not* be
degenerate; what we measured is correct.

## 4. Why $s_3$ is special and $s_{1,2}$ are not (one sentence)
$\hat n$ (normal) is globally defined, so $J\!\cdot\!\hat n$ is a scalar; the tangent plane has no global
frame, so $J_\theta,J_\phi$ are genuine vector-field components and require vector harmonics.

## 5. Fix, *if* a per-m-degenerate spatial tower is desired
Project the tangent current onto **vector spherical harmonics** instead of scalar $Y_{\ell m}\times$component:
$$\mathbf{Y}^{E}_{\ell m}=\frac{1}{\sqrt{\ell(\ell+1)}}\nabla_{S^2}Y_{\ell m},\qquad
  \mathbf{Y}^{B}_{\ell m}=\hat n\times\mathbf{Y}^{E}_{\ell m},\quad \ell\ge1,$$
and form $\Sigma^{E}_{\ell m}=\sum_n A_n\,\mathbf{Y}^{E}_{\ell m}(\hat n)\!\cdot\!\mathbf{J}_{\rm tan}(n)$
(and $E\to B$). These $E$/$B$ towers transform as proper $D^\ell$ -> icosahedral irreps -> per-m degenerate
for $\ell\le2$, splitting for $\ell\ge3$ just like $s_3$. The gradient $\nabla_{S^2}Y_{\ell m}$ in the
$(\hat\theta,\hat\phi)$ frame:
$$\nabla_{S^2}Y_{\ell m}=\partial_\theta Y_{\ell m}\,\hat\theta+\frac{1}{\sin\theta}\partial_\phi Y_{\ell m}\,\hat\phi.$$
(Conserved-current note: a conserved $J_{\rm tan}$ is divergence-free on each slice, so it is purely the
$B$/transverse part; the $E$ tower would then be the small/contact piece. Worth checking which the data wants.)

## 6. Status / recommendation
- Nothing to fix in the existing tower for correctness. For analysis, the clean rotation-invariant /
  degenerate quantities as currently coded are: $s_3$ (scalar) for $\ell\le2$, and any m-summed $g_\ell$.
- The $sp$ per-m panels are a frame/coordinate-dependent diagnostic, not irreps -- read them as such.
- Decision for NM: do we want a proper **vector-spherical-harmonic** spatial tower (Sec. 5)? That is the
  only way to get a per-m-degenerate $sp$. It is a new estimator (conn + disc), CPU-cheap to validate
  deterministically free first.

## 7. The $\ell\otimes1$ VSH expansion of the scalar-projected tangent current
$Y_{\ell m}\hat\theta$ and $Y_{\ell m}\hat\phi$ are themselves tangent **vector fields**. Adding the orbital
$\ell$ of $Y_{\ell m}$ to the spin-1 of the frame vector,
$$\ell\otimes 1=(\ell-1)\,\oplus\,\ell\,\oplus\,(\ell+1),$$
so the scalar-$Y_{\ell m}$-projected frame components expand in the **transverse** vector spherical
harmonics of total angular momentum $J=\ell-1,\ell,\ell+1$ (a tangent field has no radial part, and no
$J{=}0$ -- hairy-ball):
$$\Psi_{JM}=\nabla_{S^2}Y_{JM}\ (\text{electric/grad}),\qquad \Phi_{JM}=\hat n\times\nabla_{S^2}Y_{JM}\ (\text{magnetic/curl}),$$
with $\Phi$ contributing at $J=\ell$ (orbital $\ell$) and $\Psi$ at $J=\ell\pm1$:
$$\Sigma^{\theta}_{\ell m},\ \Sigma^{\phi}_{\ell m}=c^{B}_{\ell}\,a^{B}_{\ell m}+c^{E}_{\ell-1}\,a^{E}_{\ell-1,m}+c^{E}_{\ell+1}\,a^{E}_{\ell+1,m},$$
fixed Clebsch--Gordan $c$'s, total $m$ conserved, $a^{E/B}_{JM}=$ VSH amplitudes of $\mathbf J_{\rm tan}$.
Hence
$$sp_{\ell m}=\langle\Sigma^{\theta}\Sigma^{\theta}\rangle+\langle\Sigma^{\phi}\Sigma^{\phi}\rangle=\sum_{J,J'\in\{\ell-1,\ell,\ell+1\}} w^{(m)}_{JJ'}\,C_{JJ'}(t),$$
an **$m$-weighted (CG) superposition of the VSH two-points** $C_{JJ'}=\langle a_{Jm}a_{J'm}\rangle$ in the
$J=\ell-1,\ell,\ell+1$ channels. Each $J$ is a distinct irrep with its **own energy**, and the weights
$w^{(m)}$ are **$m$-dependent** -> per-m non-degeneracy and several effective masses (in the IR the lightest
$J$ dominates with a different weight per $m$ -> the growing splitting we measured). The scalar channel
$s_3=J\!\cdot\!\hat n$ has no such mixing (pure orbital $\ell$, $\ell\otimes0$), which is why it is clean.

Parity: every piece carries $P=(-1)^{\ell}\,(\text{from }Y_{\ell m})\times(-1)\,(\text{polar current})=(-1)^{\ell+1}$;
consistent with $\Psi$ being $(-1)^{J}$ and $\Phi$ being $(-1)^{J+1}$, all equal to $(-1)^{\ell+1}$ here. (The
**bare** local current used here is not conserved, so both $\Psi$ and $\Phi$ are present; a conserved/transverse
spatial current would be the $\Phi$ part only.)

## 8. Icosahedral ($I_h$) irrep content of the mixed $J$ channels
Each $J$ that enters branches under the proper icosahedral group $I$ (irreps $A_1, T_1, T_2, G, H$ of
dim $1,3,3,4,5$) as:

| $J$ | $I$ content | reducible (mixed)? |
| --- | --- | --- |
| 0 | $A$ | - |
| 1 | $T_1$ | no |
| 2 | $H$ | no |
| 3 | $T_2\oplus G$ | **yes** |
| 4 | $G\oplus H$ | **yes** |
| 5 | $T_1\oplus T_2\oplus H$ | **yes** |
| 6 | $A\oplus T_1\oplus G\oplus H$ | **yes** |

For $I_h$ attach the **common** parity $(-1)^{\ell+1}$ (subscript $g$ if $\ell$ odd, $u$ if $\ell$ even) to
every piece of a given $\ell$-tower -- all three $J$'s share it, so $I_h$ parity does **not** separate them;
only the rotational ($T_1,T_2,G,H$) labels do. So the scalar-projected tangent ($sp$) tower at degree $\ell$
is reducible on two counts: (a) the $J$-mixing $J=\ell{-}1,\ell,\ell{+}1$ (different energies), and (b) the
internal icosahedral splitting of any $J\ge3$. Concretely:
- $\ell=1$ ($g$): $J{=}1\,(T_{1g})\oplus J{=}2\,(H_g)$  [$J{=}0$ absent] -> 2 energy levels mixed -> not degenerate.
- $\ell=2$ ($u$): $T_{1u}\oplus H_u\oplus (T_{2u}\oplus G_u)$  ($J{=}1,2,3$; the $J{=}3$ piece mixed).
- $\ell=3$ ($g$): $H_g\oplus(T_{2g}\oplus G_g)\oplus(G_g\oplus H_g)$  ($J{=}2,3,4$; $J{=}3,4$ mixed).

Contrast the scalar tower $s_3$: a **single** $J=\ell$ channel per $\ell$, so it is fully degenerate while
$D^\ell$ is icos-irreducible ($\ell\le2$) and splits only from $\ell=3$ via $T_2\oplus G$ (Sec. 2) -- exactly
the data. To diagonalize the spatial tower into clean (per-m degenerate, $\ell\le2$) irreps, project the
tangent current directly onto $\Phi_{JM}$ and $\Psi_{JM}$ (single $J$ each), instead of scalar
$Y_{\ell m}\times$(frame component).

## 9. Continuum/CFT check: degeneracy is exact for the conserved current; bare-local is the wrong $sp$ operator
The continuum CFT formulas (qed3int_v2-14.pdf) settle whether per-m degeneracy "should" hold.

**The CFT tower is exactly m-degenerate.** Eq. (4.34/4.35) for the temporally projected correlator:
$$G^t_{\ell_1 m_1;\ell_2 m_2}(t)=C_j\big[-\tfrac13\,\delta_{\ell_1 1}\delta_{\ell_2 1}\delta_{m_1 m_2}e^{-2t}-\tfrac45\,\delta_{\ell_1 2}\delta_{\ell_2 2}\delta_{m_1 m_2}e^{-3t}+O(e^{-4t})\big].$$
The coefficients $(-\tfrac13,-\tfrac45,\dots)$ are **m-independent** and diagonal in $m$. The PDF states the
integer-spaced tower and these relative coefficients are "robust predictions from conformal theory **for the
conserved currents**." So in the continuum conformal limit all $m$ in a given $(\ell,\text{channel})$ are
exactly degenerate -- for $sp$ as well as $tp$. (Note the $\ell\ge3$ icosahedral splitting of Sec. 2/8 is a
finite-lattice artifact and is *absent* in the continuum: full $SO(3)$, every $D^\ell$ irreducible.)

**The PDF's spatial object is covariant, and uses the conserved current.** Eq. (4.27):
$$G^s(t;\hat n_1,\hat n_2)=\frac{\delta^{ab}-e^a_3 e^b_3}{D-1}\,f^{ab}_{\rm cyl}(t;\hat n_1,\hat n_2),$$
i.e. the **tangent-projected bitensor**, whose estimator Eq. (4.29) is built from the **conserved (Noether)
link current** $J^{nn'}$ (with the $\kappa$ weights), not the bare on-site $\sigma_a$. That covariant,
conserved object is what produces the degenerate tower.

**Why our cheap tower disagrees for $sp$ (and it is not a contradiction):** the local tower uses the **bare
local $\sigma_a$** projected on **scalar** $Y_{\ell m}$.
- $tp$ channel: $s_3=J\!\cdot\!\hat n$ is a genuine scalar (Sec. 4), so $s_3$ reproduces $G^t$ and is
  degenerate for $\ell\le2$ -- confirmed on the lattice.
- $sp$ channel: scalar-$Y_{\ell m}\times$(tangent frame component) is **not** the covariant $G^s$; it is the
  $J=\ell{-}1,\ell,\ell{+}1$ superposition of Sec. 7 -- a *different observable*. So it is not m-degenerate,
  and the continuum limit alone does not fix it: it is the **wrong operator for $sp$**, not a discretization error.

**Practical conclusion.**
- The local-current tower is a faithful, cheap proxy **only for $tp$ ($s_3$)**.
- For **$sp$**, use the covariant object: the **conserved link current** $G^s$ (Eq. 4.27/4.29, the expensive
  pipeline) or an equivalent **vector-spherical-harmonic** projection (Sec. 5/7). Both reproduce the
  degenerate CFT tower; the bare-local $\sigma_\theta,\sigma_\phi$ tower does not.
- On top of this, the icosahedral $\ell\ge3$ splitting is a separate $O(a)$ lattice artifact that vanishes
  in the continuum.

Refs: continuum CFT towers Eq. (4.27), (4.34/4.35), estimators (4.29)/(4.36) in qed3int_v2-14.pdf.
Vector spherical harmonics -- Barrera, Estevez, Giraldo, Eur. J. Phys. 6 (1985) 287; Thorne,
Rev. Mod. Phys. 52 (1980) 299 (electric/magnetic multipole + parity). Icosahedral $D^J$ branching and
character tables: Altmann & Herzig, *Point-Group Theory Tables*; or any $I/I_h$ correlation table.
