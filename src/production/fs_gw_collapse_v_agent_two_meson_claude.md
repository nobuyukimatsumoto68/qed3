# The GW collapse of the FS furnishing, and why $\sigma_\text{PS}=\sigma_\text{FS}$

Author: "Fin: Two-meson" session. Companion to the "{1,1,1,1}" agent's `fs_furnishing_derivation_claude.md`
(same result, independent check). Purpose: make the Ginsparg-Wilson (GW) identity and the $\sigma_\text{PS}=\sigma_\text{FS}$
statement transparent, and show exactly where the stored `peram/tau_gw` goes wrong.

## 1. The two scalar densities

From the distillation operator definitions (qed3int_v3-4.pdf Eq 5.1-5.3), a scalar density is
$$
\sigma = \eta^\dagger\, S\, \xi \;+\; \xi^\dagger\, \tilde S\, \eta ,
$$
with $S=1$ always and
$$
\tilde S = 1 \quad(\text{PS}), \qquad \tilde S = -(1-D_\text{ov}^\dagger) \quad(\text{FS}).
$$
So $\sigma_\text{PS}=\eta^\dagger\xi+\xi^\dagger\eta$ (plain), and $\sigma_\text{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_\text{ov}^\dagger)\eta$
(the second term carries the "furnishing" operator $(1-D_\text{ov}^\dagger)$). The two look different; the claim is
that inside any correlator they give the SAME thing.

## 2. The toolbox: two exact identities of the massless overlap

- **Ginsparg-Wilson (IV.17 of qed3_v2-6.pdf):**
$$
D_\text{ov}^{-\dagger} \;=\; 1 - D_\text{ov}^{-1}, \qquad\text{equivalently}\qquad D_\text{ov}^{-1}+D_\text{ov}^{-\dagger}=1 .
$$
This is the massless GW relation written for the inverses. It says the BACKWARD (adjoint) propagator is
$1$ minus the FORWARD propagator.
> **CORRECTION (2026-09-17, NM):** this system does NOT have $\sigma_3$- (gamma-) hermiticity -- in 3D / 2-component
> there is no chirality and $D_\text{ov}^{-\dagger}\neq\sigma_3 D_\text{ov}^{-1}\sigma_3$. An earlier draft here (and
> the `{1,1,1,1}` thread's memory) claimed the suppression $\langle\sigma_\text{PS}\,\sigma_\text{PS}^2\rangle=0$ was
> "protected by $\sigma_3$-hermiticity ($\mathrm{Tr}[\Gamma_\text{even}GGG]=0$)". **That mechanism is WRONG and is
> retracted.** What is empirically true (directly verified): $\langle\sigma^2|m_{PS}\rangle\approx0$ at the COMPLETE
> distillation basis (L1 Nv=24=2$\cdot$12, L2 Nv=84, and free), and it LEAKS once the basis is truncated (L2 Nv=24 of
> 84) -- see the free-case test in `project_sigma_sigma_f2_mixing.md`. The correct mechanism is not yet derived; do
> NOT attribute it to $\sigma_3$-hermiticity. NOTE: everything below (the collapse, PS$=$FS, the `tau_gw` bug) rests
> ONLY on the Ginsparg-Wilson relation above, NOT on $\sigma_3$-hermiticity, so it stands.

Write the forward propagator as $G\equiv D_\text{ov}^{-1}$ (this is exactly the plain perambulator leg
$\tau=V^\dagger G V$). Its adjoint is $G^\dagger=D_\text{ov}^{-\dagger}$.

## 3. The collapse

The FS furnishing $\tilde S=-(1-D_\text{ov}^\dagger)$ sits on the daggered ($\xi^\dagger\ldots\eta$) side of the
vertex, so in the Wick contraction it multiplies the leg's **adjoint** (backward) propagator $G^\dagger=D_\text{ov}^{-\dagger}$
-- NOT the forward $G$. Now compute the furnished leg, using only $D_\text{ov}^\dagger D_\text{ov}^{-\dagger}=1$
and GW:
$$
\tilde S\,G^\dagger
= -(1-D_\text{ov}^\dagger)\,D_\text{ov}^{-\dagger}
= -\big(D_\text{ov}^{-\dagger} - \underbrace{D_\text{ov}^\dagger D_\text{ov}^{-\dagger}}_{=\,1}\big)
= -\big(D_\text{ov}^{-\dagger}-1\big)
= 1 - D_\text{ov}^{-\dagger}
\;\overset{\text{GW}}{=}\; D_\text{ov}^{-1} = G .
$$
So the furnished FS leg **collapses to the plain forward propagator** $G=\tau$. The $-(1-D_\text{ov}^\dagger)$
exactly converts the backward propagator back into the forward one. Since the PS leg is already $G$, the FS and
PS densities produce IDENTICAL contractions:
$$
\boxed{\;\sigma_\text{FS}\ \cong\ \sigma_\text{PS}\ \text{in every correlator}\;\Longrightarrow\; \text{PS}=\text{FS}.\;}
$$

### Why, physically
The FS density is built from the "backward/adjoint" leg dressed by $(1-D_\text{ov}^\dagger)$. At $m=0$ the GW
relation makes "backward $=1-$ forward", and the furnishing is precisely the operator that undoes the "$1-$".
So FS is not a genuinely new operator -- it is the same forward-propagator object written in a backward-leg
disguise. This is the operator-level statement behind the measured $\text{PS}=\text{FS}$ (scalar 2pt) and
$\text{PP}=\text{FF}$ (the $\sigma^2$ four-point, bit-identical). The flavor factor
$\prod_\text{loops}(1+(-1)^{n_\text{FS}})$ used in the `flavorgeom` cache is just the combinatorial bookkeeping
of this same collapse (each FS vertex either reproduces the PS one or contributes a sign that cancels on the
tadpole loops, which vanish under the contact subtraction).

## 4. Where `peram/tau_gw` goes wrong (the bug)

The stored second perambulator is built (`distill_peram_claude.cu:463-473`) as
$$
\texttt{peram/tau\_gw} \;=\; V^\dagger\,(1-D_\text{ov}^\dagger)\,D_\text{ov}^{-1}\,V ,
$$
i.e. it applies $(1-D_\text{ov}^\dagger)$ to the **forward** inverse $D_\text{ov}^{-1}$ (a forward solve, then one
$(1-D_\text{ov}^\dagger)$ mat-vec), NOT to the adjoint inverse. That is the wrong leg. Using
$D_\text{ov}^\dagger D_\text{ov}^{-1}=D_\text{ov}^\dagger-1$ (multiply GW $D_\text{ov}^{-1}+D_\text{ov}^{-\dagger}=1$
on the left by $D_\text{ov}^\dagger$):
$$
(1-D_\text{ov}^\dagger)D_\text{ov}^{-1}
= D_\text{ov}^{-1} - (D_\text{ov}^\dagger-1)
= 1 + D_\text{ov}^{-1} - D_\text{ov}^\dagger .
$$
This carries the **bare** $D_\text{ov}^\dagger$ (the Dirac adjoint OPERATOR, an $O(1)$ object), whereas the
collapse needs the adjoint INVERSE $D_\text{ov}^{-\dagger}$. So $\tau_{gw}\neq\tau$; the difference
$\tau_{gw}-\tau = -(D_\text{ov}^\dagger-1)-\ldots$ is exactly the spurious bare-$D_\text{ov}^\dagger$ piece that
faked an "$O(1)$ FS coupling" in the sigma^2 triangle/four-point.

### The correct object is free
The adjoint inverse is available with NO stored `tau_gw` and NO new solve, straight from GW:
$$
D_\text{ov}^{-\dagger} \;=\; 1 - D_\text{ov}^{-1} \;\Longrightarrow\; \bar\tau \equiv V^\dagger D_\text{ov}^{-\dagger}V = \delta - \tau = 1 - \tau .
$$
And since $\tilde S\,\bar\tau$ collapses to $\tau$ (Sec 3), the FS leg is simply $\tau$. So: **use $\tau$ for the
FS leg** (or, if you want the adjoint kept explicit, $\delta-\tau=1-\tau$), and **never `tau_gw`**.

## 5. Checks

- Algebra: two independent derivations (this note; `fs_furnishing_derivation_claude.md`) agree.
- Suppression $\langle\sigma_\text{PS}(t)\,\sigma_\text{PS}^2(0)\rangle\approx0$: empirically holds at the COMPLETE
  basis (mechanism NOT $\sigma_3$-hermiticity -- see the CORRECTION note in Sec 2; the system lacks it). It LEAKS
  under basis truncation. With the collapse, FS behaves as PS, so the same completeness-dependence applies.
- Numerical (free L1): the sigma^2 triangle with the $\tau$ (collapsed) leg gives $R\sim10^{-6}=0$; with
  $-\tau_{gw}$ it gave a spurious $-0.0046$. Production `four_point` FS.FS, after replacing $-\tau_{gw}\to\tau$,
  is now $==$ PS.PS exactly ($C_S$ effmass $0.3891$), where before it was $\sim0.52\times$PS.PS.

## 6. One-line takeaways

- GW: $D_\text{ov}^{-1}+D_\text{ov}^{-\dagger}=1$. Backward propagator $=1-$ forward.
- FS furnishing on the (correct) adjoint leg collapses: $-(1-D_\text{ov}^\dagger)D_\text{ov}^{-\dagger}=D_\text{ov}^{-1}=\tau$. Hence PS$=$FS.
- `tau_gw` applied $(1-D_\text{ov}^\dagger)$ to the FORWARD inverse instead $\Rightarrow$ carries bare $D_\text{ov}^\dagger$ $\Rightarrow$ artifact.
- Fix: FS leg $=\tau$ (or $\delta-\tau$), never `tau_gw`.
