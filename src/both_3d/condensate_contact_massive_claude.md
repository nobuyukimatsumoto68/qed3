# Isolating the contact term for the massive overlap inverse ($m_P$, $m_F$)

Goal: reproduce the Luscher-algebra contact-term isolation of `general_overlap_v3.pdf` Eqs. (1.7)/(1.15)
for the massive operators of `qed3int_v2(7)/main.tex` (= qed3int_v2-14.pdf) Sec. "Pseudofermion ... soft
mass" (Eqs. around 975-997). NOTE: use the current draft `qed3int_v2(7)/`, NOT the stale root `main.tex`.

Lattice units throughout (the note's convention). Parametrization (note 1.1-1.2, 1.8):
$$
D=1+V,\qquad V=\frac{i\gamma\cdot\tilde p+b}{w},\qquad \tilde p^2+b^2=w^2\ (\Leftrightarrow V\ \text{unitary}).
$$
Continuum/small-$p$: $w\to1$, $b\to-1$, $\tilde p\to p$, $w+b\to\tfrac12 p^2$ (note 1.14).

## 0. Massless recap (note 1.7 / 1.15)
$$
D^{-1}=\tfrac12\Big(1+\frac{-i\gamma\cdot\tilde p}{w+b}\Big)
=\underbrace{\tfrac12}_{\text{contact}}+\underbrace{\tfrac12\,\frac{-i\gamma\cdot\tilde p}{w+b}}_{\text{propagating}}
\ \xrightarrow{a\to0}\ \tfrac12+\frac{-i\gamma\cdot p}{p^2}\,(1+O(a^2)).
$$
The contact $\tfrac12\,\mathbb 1$ is the momentum-INDEPENDENT (i.e. $\delta_{xy}$) part. It is fixed by GW
(note 1.6): writing $D^{-1}=c\,\mathbb 1+(\gamma_5\text{-odd})$, the relation $\{D^{-1},\gamma_5\}=\gamma_5$ gives
$2c=1$, so $c=\tfrac12$ exactly. (Matches your free data: $\operatorname{Re}\operatorname{tr}[A\,D_{ov}^{-1}]/\operatorname{tr}[A\mathbb 1]=\tfrac12$.)
Note: `general_overlap_v3.pdf` (1.15) prints "$+1$"; the consistent value is "$+\tfrac12$" (cf. 1.7), unless
that note normalizes $D\to 2D$.

## 1. Master lemma: $(D_{ov}+m)^{-1}$ for any $m\in\mathbb C$
$D+m=\dfrac{(1+m)w+b+i\gamma\cdot\tilde p}{w}$. Rationalize using $(\,S+i\gamma\tilde p)(S-i\gamma\tilde p)=S^2+\tilde p^2$
with $S=(1+m)w+b$, and $b^2+\tilde p^2=w^2$:
$$
(D_{ov}+m)^{-1}=\frac{[(1+m)w+b]-i\gamma\cdot\tilde p}{\Delta_m},\qquad
\boxed{\ \Delta_m\equiv[(1+m)^2+1]\,w+2(1+m)\,b=m^2 w+2(1+m)(w+b).\ }
$$
Now split off the SAME constant $\tfrac12$. Since $[(1+m)w+b]-\tfrac12\Delta_m=-\tfrac12\,m(mw+2b)$,
$$
\boxed{\ (D_{ov}+m)^{-1}=\tfrac12\;-\;\frac{m\,(mw+2b)}{2\,\Delta_m}\;+\;\frac{-i\gamma\cdot\tilde p}{\Delta_m}.\ }
$$
- Term 1: the contact $\tfrac12\,\mathbb 1$ -- momentum-independent and **mass-independent**, identical to massless.
- Term 2: a momentum-DEPENDENT scalar ($\propto\mathbb 1$ in spin). It is NOT a contact term; it is the genuine
  scalar part of the massive propagator. It is $O(m)$ and vanishes as $m\to0$.
- Term 3: the propagating $\gamma$-odd part.

Check $m=0$: $\Delta_0=2(w+b)$, terms 2,3 $\to$ note (1.7). 

### Continuum form (massive analog of 1.15)
$\Delta_m\to m^2+(1+m)p^2$, $\;mw+2b\to m-2$, so
$$
(D_{ov}+m)^{-1}\ \xrightarrow{a\to0}\ \tfrac12+\frac{m-i\gamma\cdot p}{p^2+m^2}\,(1+O(a^2))
=\tfrac12+\frac{1}{\,i\gamma\cdot p+m\,}(1+O(a^2)).
$$
i.e. contact $\tfrac12$ plus the standard Euclidean propagator. (Scheme caveat: with $D_{ov}+m$ the residue
carries an extra $O(m^2)$ piece, $\frac{m-i\gamma p}{p^2+m^2}-\frac{m^2/2}{p^2+m^2}$; the properly normalized
massive overlap $(1-\tfrac m2)D_{ov}+m$ removes it. Your code uses $D_{ov}+m$, so keep the exact boxed form.)

## 2. Parity-symmetric case ($\sigma_{PS}$, real $m=m_F$)
(`qed3int_v2(7)/main.tex:975,1013`: parity-symmetric / flavor-breaking term $m_F\sigma_{PS}$, $m_F\in\mathbb R$.)
Both legs are $(D_{ov}+m_F)^{\pm}$ (backward $=$ dagger, since $m_F\in\mathbb R$). Use the lemma directly:
$$
(D_{ov}+m_F)^{-1}=\tfrac12-\frac{m_F(m_F w+2b)}{2\Delta_{m_F}}+\frac{-i\gamma\cdot\tilde p}{\Delta_{m_F}} .
$$
Contact $=\tfrac12$, real. Scalar part real -> the parity-symmetric condensate $\operatorname{tr}$ of term 2.

## 3. Flavor-symmetric case ($\sigma_{FS}$, imaginary $m=m_P$)
(`qed3int_v2(7)/main.tex:986,1013`: flavor-symmetric / parity-breaking term $m_P\sigma_{FS}$, $m_P\in i\mathbb R$.)
**Forward leg** $(D_{ov}+m_P)^{-1}$: the lemma verbatim with $m=m_P$ (imaginary). Contact still $\tfrac12$;
scalar part term 2 is now imaginary ($m_P^2<0$).

**Backward leg** $B\equiv(1+m_P)D_{ov}^\dagger-m_P$ (`:992-993`). Factor $B=(1+m_P)\,[\,D_{ov}^\dagger-\tilde m_P\,]$,
$\;\tilde m_P\equiv\dfrac{m_P}{1+m_P}$. Apply the lemma to $D^\dagger=1+V^\dagger$ (i.e. $\tilde p\to-\tilde p$,
mass $-\tilde m_P$):
$$
B^{-1}=\frac{1}{1+m_P}\Big[\ \tfrac12\;-\;\frac{\tilde m_P(\tilde m_P w-2b)}{2\,\tilde\Delta}\;+\;\frac{+i\gamma\cdot\tilde p}{\tilde\Delta}\ \Big],
\qquad \tilde\Delta\equiv \tilde m_P^2 w+2(1-\tilde m_P)(w+b).
$$
So the **flavor backward leg has a RESCALED contact** $\dfrac{1}{2(1+m_P)}\,\mathbb 1$ (not $\tfrac12$), plus
the prefactor $\frac{1}{1+m_P}$ on the rest. This $\frac{1}{1+m_P}$ is the same anomalous factor that appears
in the reweighting $R$ (`qed3int_v2(7)/main.tex:1016-1036`). In the continuum $\tilde m_P\to m_P$, contact $\to\tfrac12$.
(The memory's $\tilde D_{m_P}=D_{ov}+\tfrac{m_P}{1-m_P}$ is the conjugated/forward form of this same $B$.)

## 4. Convention (current draft = qed3int_v2-14 = qed3int_v2(7)/main.tex)
The subscript labels the symmetry that the mass term BREAKS (`:975,986,1013`):
- $m_F\in\mathbb R$ multiplies $\sigma_{PS}$ -> parity-symmetric, **flavor-breaking** case; backward leg $=$ dagger.
- $m_P\in i\mathbb R$ multiplies $\sigma_{FS}$ -> flavor-symmetric, **parity-breaking** case; backward leg $(1+m_P)D^\dagger-m_P$.
This MATCHES the code/memory ($m_F$ real/clean, $m_P$ imaginary/$\tilde D$). (The stale root `main.tex` used the
opposite pairing $m_P\sigma_P$/$m_F\sigma_F$ AND had a stray $m_P$-for-$m_F$ slip in the flavor forward leg --
disregard that file.) Contact summary:
- real $m_F$, $\sigma_{PS}$: contact $\tfrac12$ both legs.
- imag $m_P$, $\sigma_{FS}$: contact $\tfrac12$ (fwd), $\frac{1}{2(1+m_P)}$ (bwd).

## 5. Consequence for the condensate divergence
The $\delta_{xy}$ (momentum-independent) divergence in $\langle\bar\psi\psi\rangle$ is exactly the contact:
- $\sigma_{PS}=-2\operatorname{Re}\operatorname{tr}D_m^{-1}$: subtract $-2\cdot\tfrac12\,\mathbb 1=-\mathbb 1$
  (i.e. measure $-2\operatorname{Re}\operatorname{tr}[A(D_m^{-1}-\tfrac12)]$), removing the $m$-independent contact.
- The GW-improved operator $\bar\psi(1-\tfrac12 D_{ov})\psi$ does this automatically (cf.
  `condensate_subtraction_claude.md`). NB the flavor combination $(1-D_{ov})$ in $\sigma_{FS}$ uses the FULL
  $V^\dagger$, not $\tfrac12$, so it does NOT cancel the contact (gives $-\operatorname{tr}\mathbb 1$ at $m=0$);
  the $\tfrac12$ subtraction is still needed.
- Caveat (subleading): term 2 summed over the BZ is itself linearly divergent $\sim m\Lambda$ in 3D, so the
  $\tfrac12$-subtraction kills only the leading ($m$-independent) divergence; the residual $O(m)$ divergence
  needs the free-field subtraction or the $(1-\tfrac m2)$-improved operator, as before.

## 6. Interacting generalization (gauge field on)
The free derivation diagonalized in MOMENTUM. The only things it actually used are (i) $V\equiv D_{ov}-1$
UNITARY and (ii) $\gamma_5$-Hermiticity -- both hold for ANY gauge background (they ARE the overlap/GW
construction, exact up to the Zolotarev sign-approx error). So the contact survives verbatim.

### (a) Massless contact is exact and config-independent
From $V$ unitary alone ($V^\dagger=V^{-1}$):
$$
D_{ov}^{-1}+D_{ov}^{-\dagger}=(1+V)^{-1}+(1+V^{-1})^{-1}=(1+V)^{-1}+(1+V)^{-1}V=(1+V)^{-1}(1+V)=\mathbb 1 .
$$
Hence the Hermitian part of $D_{ov}^{-1}$ is $\tfrac12\mathbb 1$ EXACTLY, on every configuration -- no momentum
space needed. So $\operatorname{Re}\operatorname{tr}[A\,D_{ov}^{-1}]=\tfrac12\operatorname{tr}[A]=V_{st}$ for free AND
interacting alike. The $\tfrac12$ contact you subtract is literally the same operator; the free run is not
just illustrative, it exhibits the exact interacting contact.

### (b) Momentum $\to$ overlap eigenphase
Replace the c-number $\frac{-i\gamma\cdot\tilde p}{w+b}$ by the operator $K\equiv 2D_{ov}^{-1}-1=\frac{1-V}{1+V}$
(Cayley transform of $V$; anti-Hermitian since $V$ unitary: $K^\dagger=-K$). The eigenvalues of $D_{ov}$ sit on
the GW circle $\lambda_i=1+e^{i\theta_i}$ -- the SAME circle free and interacting; the gauge field only
redistributes the phases $\theta_i$ ALONG it, never off it. Per eigenvalue, identically to the master lemma
with $w\to1,\ b\to\cos\theta_i$ (so $mw+2b\to m+2\cos\theta_i$, $\Delta_m\to D_{\theta_i}$):
$$
\operatorname{Re}\frac{1}{\lambda_i+m}=\tfrac12-\frac{m(m+2\cos\theta_i)}{2\,D_{\theta_i}},
\qquad D_{\theta_i}=(1+m)^2+2(1+m)\cos\theta_i+1 .
$$
Therefore, EXACTLY (any config), with $\{\theta_i\}$ the eigenphases of $V=D_{ov}-1$:
$$
\boxed{\ \operatorname{Re}\operatorname{tr}D_m^{-1}=\frac{N}{2}\ -\ \sum_i\frac{m(m+2\cos\theta_i)}{2\,D_{\theta_i}}\ }
$$
contact $\tfrac N2$ (config-independent) $+$ dynamical, gauge-dependent term ($O(m)$, vanishes as $m\to0$).
This is the interacting analog of the lemma's "term 1 $+$ term 2"; the free case is just the special
spectrum $\theta_i=\theta(p_i)$.

### (c) What does NOT transfer, and why free-subtraction is justified
The explicit closed form ($w,b,\tilde p$ as functions of one momentum) does not, because $V$ is not
momentum-diagonal -- only the EIGENPHASE formula above is general. Two consequences:
- The residual $O(m)$ divergence (the UV tail of $\sum_i$) comes from eigenphases on the DENSE/UV part of the
  circle; those UV modes are gauge-insensitive (asymptotic freedom), so the divergent tail is the same free
  and interacting -> the free-field subtraction removes it EXACTLY up to the IR. This is the rigorous
  justification of the free subtraction in `condensate_subtraction_claude.md`.
- The IR eigenphases ($\lambda_i\to0$, $\theta_i\to\pi$) carry the dynamics: as $m\to0$ the dynamical sum is
  governed by the spectral density at the origin -> Banks-Casher / the order parameter (the Karthik-Narayanan
  spectral route, `condensate_subtraction_claude.md` Sec. 5). Flavor case ($m_P$ imag, $\sigma_{FS}$): the
  same per-eigenphase formula with $m\to m_P$ and the backward leg's $\frac{1}{1+m_P}$ prefactor.

## 7. Is the $(1-\tfrac12 D)$ subtraction "just contact removal", or deeper?
Deeper. Contact removal is the visible COROLLARY; the content is Luscher's EXACT lattice chiral symmetry.

**(i) Dictated by an exact symmetry.** GW fermions have an exact symmetry of the lattice action
(Luscher, hep-lat/9802011), and it is the one in this paper (`main.tex` eq:chiral_symmetry):
$$
\delta\psi=\hat\gamma_5\,\psi,\quad \delta\bar\psi=\bar\psi\,\gamma_5,\quad \hat\gamma_5\equiv\gamma_5(1-D_{ov}).
$$
The rotation uses $\hat\gamma_5$ on $\psi$ but $\gamma_5$ on $\bar\psi$. The symmetric generator is
$\tfrac12(\gamma_5+\hat\gamma_5)=\gamma_5(1-\tfrac12 D_{ov})$, so $\bar\psi(1-\tfrac12 D_{ov})\psi$ is the scalar
density on which the exact symmetry closes as an honest mass term. The factor $\tfrac12$ is FIXED by the
symmetry, not chosen to cancel a number.

**(ii) Propagator form of the same fact.** The GW relation is a defect in the continuum chiral algebra,
$\{\gamma_5,D_{ov}^{-1}\}=\gamma_5$ (note Eq. 1.6), while the continuum massless propagator obeys
$\{\gamma_5,S\}=0$. The subtraction restores it exactly:
$$
\{\gamma_5,\;D_{ov}^{-1}-\tfrac12\}=\gamma_5-\tfrac12\cdot2\gamma_5=0 .
$$
So the $\tfrac12$ contact IS the GW defect; subtracting it restores the continuum chiral anticommutation of
the propagator, identically on every config (cf. Sec. 6a).

**(iii) The payoff bare contact-removal would not give.** Protected by the exact symmetry,
$\bar\psi(1-\tfrac12 D_{ov})\psi$ does NOT mix additively with $\mathbb 1$: it is MULTIPLICATIVELY
renormalizable ($Z_S$ only) and obeys the correct chiral Ward identities, so it has a genuine continuum limit
as an order parameter. The naive $\bar\psi\psi$ mixes with $\sim a^{-2}\mathbb 1$.

**(iv) 3D caveat -- two distinct $(1-D)$ objects (do not conflate).**
- $(1-\tfrac12 D_{ov})$: the Niedermayer/Hasenfratz/Luscher contact subtraction (factor $\tfrac12$) -> chiral
  covariance / contact removal. NEEDED for the order parameter in BOTH mass cases.
- full $(1-D_{ov})=-V^\dagger$ inside $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta$: the
  flavor-multiplet projector ($\hat P_\pm$, eq:projectors), making $\sigma_{FS}=\sum_i\bar\psi_i\psi_i$ the
  exact-FLAVOR-symmetric singlet. This is NOT a contact subtraction -- it gives $-\operatorname{tr}\mathbb 1$ at
  $m=0$, not $0$.

Same GW origin, different jobs: full $(1-D)$ enforces exact lattice FLAVOR symmetry; the $\tfrac12$ enforces
chiral covariance / contact removal. Refs: Luscher hep-lat/9802011; Hasenfratz-Laliena-Niedermayer
hep-lat/9801021; see `condensate_subtraction_claude.md` Sec. 5.

## 8. Data check: do these formulas explain the $m_F=0.1$ shifts from $-2$? (YES)
Free $L=1$, dense ref (`jj_validate_mF_claude.ipynb` cell 27), real $m_F=0.1$. Everything is REAL here
(`etadag_xi` imag $\sim10^{-18}$). Write $T\equiv\operatorname{tr}[A\,D_{m_F}^{-1}]=-$`etadag_xi`,
$N\equiv\operatorname{tr}[A\mathbb 1]=2V_{st}$. Data:
$$
T/V_{st}=0.96192,\qquad \sigma_{PS}/V_{st}=-1.9238,\qquad \sigma_{FS}/V_{st}=-1.9038 .
$$

### Both condensates are ONE number $T$ (not two independent observables)
- $\sigma_{PS}=-2\,\operatorname{Re}T$ (definition).
- $\sigma_{FS}=-T+\overline{(1+m_F)T-N}$ via the EXACT identity $(1-D_{ov})D_m^{-1}=(1+m)D_m^{-1}-\mathbb 1$
  (just $D_m=D_{ov}+m$, no GW). For real $m_F$: $\operatorname{Re}\sigma_{FS}=m_F\,\operatorname{Re}T-N$.

So the only dynamical input is $T$; PS and FS are two linear maps of it.

### The shift is term 2 (mass-dependent), NOT the contact
Define $\delta T\equiv \operatorname{Re}T-V_{st}$ (the term-2 / dynamical movement of the trace off its
contact value $V_{st}$). Here $\delta T/V_{st}=0.96192-1=-0.03808$. Then:
$$
\frac{\sigma_{PS}}{V_{st}}=-2-2\,\frac{\delta T}{V_{st}}=-2+0.0762,
$$
$$
\frac{\sigma_{FS}}{V_{st}}=-2+\underbrace{m_F}_{\text{contact-driven}}+\underbrace{m_F\frac{\delta T}{V_{st}}}_{\text{dynamical}}
   =-2+\underbrace{0.1}_{}\underbrace{-0.0038}_{}=-2+0.0962 .
$$
Both match the data to all displayed digits. Note the asymmetry:
- the PS shift ($+0.0762$) is PURELY dynamical ($-2\,\delta T$);
- the FS shift ($+0.0962$) is MOSTLY the trivial contact-driven piece $m_F\cdot V_{st}$ ($=+0.1$), with only a
  small dynamical correction $m_F\,\delta T$ ($=-0.0038$). I.e. $\sigma_{FS}$'s shift is dominated by
  $m_F\times$(the $\tfrac12$ contact), not by physics.

### Equivalent: the FS-PS relation (eliminate $T$)
$$
\boxed{\ \sigma_{FS}=-\tfrac{m_F}{2}\,\sigma_{PS}-2V_{st}\ }\qquad(\text{real }m_F,\text{ free}).
$$
Check: $-0.05(-1.9238)-2=-1.9038$. A pure consistency relation -- no new information in $\sigma_{FS}$.

### What the formulas DO and DO NOT do
- DO: reproduce both shifted numbers and their relation, from the single $T$.
- DO NOT (by themselves): split $\delta T$ into physical condensate vs residual $O(m)$ UV artifact. By §6,
  $\delta T=-\sum_i\frac{m_F(m_F+2\cos\theta_i)}{2D_{\theta_i}}$ (area-weighting aside) mixes both; the IR
  eigenphases are the order parameter, the UV tail is the artifact removed by free-subtraction.

CAVEAT (area weighting): the eigenphase sum in §6 is for the UNWEIGHTED $\operatorname{tr}D_m^{-1}$; the data
uses $\operatorname{tr}[A\,D_m^{-1}]$. The contact part is unaffected ($\operatorname{Re}\operatorname{tr}[A\,D_{ov}^{-1}]=\tfrac12\operatorname{tr}[A]$
exactly), but a literal $\sum_i\frac{1}{\lambda_i+m}$ reproduces $\delta T$ only if $A$ commutes with $D_{ov}$
(true per-site for the contact, approximate for term 2). To check $\delta T$ from first principles, evaluate
$\operatorname{Re}\operatorname{tr}[A(D_{m_F}^{-1}-\tfrac12)]$ on the dense free $D_{ov}$ directly.

## 9. The small cross-channel shift in the OTHER condensate (and a csub slip)
Observation (`jj_validate_{mF,mP}`): after contact subtraction, the EXPECTED channel shifts $O(m)$, but the
OTHER channel still shows a shift $\sim10\times$ smaller. Two contributions:

### (a) NO genuine cross-channel order parameter -- symmetry protects it
[Correction of an earlier sloppy "leakage" framing.] The mF run action $m_F\sigma_{PS}$ is PARITY-INVARIANT
($\sigma_{PS}$ parity-even, $m_F$ real). The genuine order parameter in the OTHER channel is the PARITY-ODD
part of $\sigma_{FS}$; by EXACT lattice parity its VEV is ZERO. (Already visible at $m=0$:
$\langle\sigma_{FS}\rangle=-N\neq0$ in the parity-symmetric free theory -- so that nonzero value is the
PARITY-EVEN admixture of $\sigma_{FS}$, NOT an order parameter.) So there is no symmetry-violating leakage.
The observed $\sigma_{FS}$ shift lives entirely in the parity-EVEN component:
- DOMINANT: the MASS-DEPENDENT contact (see (c)) -- $\sigma_{FS}$'s momentum-independent part is $(m_F-2)V_{st}$,
  shifting $+m_F V_{st}$ from $-2V_{st}$.
- small: a parity-even dynamical residual $m_F\delta T=-0.0038$ ($\delta T\equiv\operatorname{Re}T-V_{st}$),
  $=-\tfrac{m_F}{2}\times(\sigma_{PS}\text{ shift})$.
Both are parity-even and $\to0$ in the continuum. (The mP run DIFFERS: there the cross channel $\sigma_{PS}$ is
protected by FLAVOR, which is ANOMALOUS on the lattice -- so a cross shift there need NOT vanish by symmetry;
cf. the parity-anomaly / Chern-Simons / $R$ discussion in the paper.)

### (b) Contact-subtraction slip in the notebook (per-RUN vs per-observable)
The rescaled backward-leg contact $\frac{1}{2(1+m)}$ (§3) is a property of the RUN (which mass is in the
ACTION), not of the observable. mP run (flavor action) = rescaled backward leg on ALL observables; mF run
(parity action) = CLEAN dagger backward leg on all observables. Correct $\delta_{xy}$-contact (csub):

| run | $\sigma_{PS}$ csub | $\sigma_{FS}$ csub |
|---|---|---|
| mF (clean bwd)     | $2$                  | $2-m_F$ |
| mP (rescaled bwd)  | $1+\tfrac{1}{1+m_P}$ | (rederive: $\tilde D_{m_P}$+$(1-D)$; NOT $2$) |

`jj_validate_*` cells use `csub = (1+1/(1+m)) if MASS=='mP'(PS) / 'mF'(FS) else 2.0`. This is correct for
mF-PS and mP-PS, but WRONG for mF-FS (uses $1+\tfrac{1}{1+m_F}$ instead of $2-m_F$; over-subtracts by
$\frac{m_F^2}{1+m_F}=+0.0091$) and for mP-FS (uses $2$, missing the $m_P$-dependence). The mF-FS residual
decomposes as $+0.0053=\underbrace{+0.0091}_{\text{csub slip}}\underbrace{-0.0038}_{\text{genuine leakage (a)}}$.
FIX: csub by RUN (backward leg) and add the $(1-D)$-induced $m$-shift for $\sigma_{FS}$. mF: PS$\to2$,
FS$\to2-m_F$. mP: PS$\to1+\tfrac{1}{1+m_P}$, FS$\to$ derive. Then the cross-channel residual is purely (a).

### (c) The contact term IS mass-dependent (the simplest correct explanation)
Momentum-independent ($\delta_{xy}$) part of each leg, closed form:
- plain $D_m^{-1}$: $\tfrac12$  (mass-INDEPENDENT).
- $(1-D_{ov})D_m^{-1}=(1+m)D_m^{-1}-1$: $\tfrac{m-1}{2}$  (mass-DEPENDENT).
- rescaled backward $[(1+m)D^\dagger-m]^{-1}$: $\tfrac{1}{2(1+m)}$  (mass-DEPENDENT).

Hence $\sigma_{PS}$ (both legs plain) has a mass-INDEPENDENT contact $-2V_{st}$ -> its WHOLE shift is dynamical
($-2\delta T$). $\sigma_{FS}$ (one leg carries the $(1-D)$ dressing) has a mass-DEPENDENT contact $(m_F-2)V_{st}$
-> its shift is DOMINATED by the contact mass-dependence $+m_F V_{st}$, plus the small $m_F\delta T$. This is
the easiest correct picture of the shifts. The exact full-order contact subtraction = subtract these
closed-form m.i. parts per channel (this is what csub should be).

CAUTION on the §5/§7 $(1-\tfrac12 D_{ov})$ operator (the other writer's suggestion): inserting it does NOT
remove the contact in full order at finite mass. $(1-\tfrac12 D_{ov})D_m^{-1}=(1+\tfrac m2)D_m^{-1}-\tfrac12$
has m.i. part $\tfrac12-\tfrac12(1-\tfrac m2)=\tfrac m4\neq0$: it kills the $m=0$ contact but leaves an $O(m)$
contact $\tfrac m4$ (for $D_m=D_{ov}+m$; normalization-dependent). Its justification is operator
renormalization / chiral covariance (§7), NOT the chiral limit -- so it is not WRONG to use with explicit
SSB-probing mass, but it is neither necessary nor a complete contact subtraction here. To understand/clean the
shifts, subtract the exact m.i. parts above; to take $m\to0$ for the SSB check, the contact is analytic in $m$
and drops out of the non-analytic discontinuity regardless.

## 10. Exact per-channel contact subtraction (all four csub) -- derivation for hand-check
Free $L=1$; the m.i. (ultralocal $\delta_{xy}$) contact rule $[\![(D_{ov}+\mu)^{-1}]\!]=\tfrac12$ is EXACT in the
free theory (master lemma term 1, momentum-independent), so every csub below is exact for the free tests.
$N\equiv\operatorname{tr}[A\mathbb 1]=2V_{st}$; contact of $\operatorname{tr}[AX]$ in units of $V_{st}$ is $2[\![X]\!]$.

### Contractions (per run)
Action $S=\eta^\dagger M_f\xi+\xi^\dagger M_b\eta$ -> forward prop $\langle\xi\eta^\dagger\rangle=M_f^{-1}$,
backward prop $\langle\eta\xi^\dagger\rangle=M_b^{-1}$. With the Grassmann reorder sign:
$$
\texttt{etadag\_xi}=-\operatorname{tr}[A\,M_f^{-1}],\quad
\texttt{xidag\_eta}=-\operatorname{tr}[A\,M_b^{-1}],\quad
\texttt{xidag\_1mDdag\_eta}=-\operatorname{tr}[A\,(1-D_{ov}^\dagger)\,M_b^{-1}],
$$
$\sigma_{PS}=$`etadag_xi`$+$`xidag_eta`, $\sigma_{FS}=$`etadag_xi`$-$`xidag_1mDdag_eta`. Per run
(`qed3int_v2(7)/main.tex:975-996`):
- mF (real $m_F$): $M_f=D_{ov}+m_F$, $M_b=D_{ov}^\dagger+m_F$.
- mP (imag $m_P$): $M_f=D_{ov}+m_P$, $M_b=(1+m_P)D_{ov}^\dagger-m_P=(1+m_P)(D_{ov}^\dagger+\nu)$,
  $\nu\equiv-\tfrac{m_P}{1+m_P}$. [Consistency: $(D_{ov}^\dagger+\nu)^{-1}=\tilde D_{m_P}^{-\dagger}$ since
  $\overline{m_P/(1-m_P)}=\nu$ for imaginary $m_P$ -- matches `phimm`.]

### Building-block contacts
$[\![(D_{ov}+\mu)^{-1}]\!]=\tfrac12$ (any $\mu$, daggered too). And from
$(1-D_{ov}^\dagger)(D_{ov}^\dagger+\mu)^{-1}=(1+\mu)(D_{ov}^\dagger+\mu)^{-1}-1$:
$[\![(1-D_{ov}^\dagger)(D_{ov}^\dagger+\mu)^{-1}]\!]=\tfrac{\mu-1}{2}$. Backward-leg contacts:
- mF: $[\![M_b^{-1}]\!]=\tfrac12$; $[\![(1-D_{ov}^\dagger)M_b^{-1}]\!]=\tfrac{m_F-1}{2}$.
- mP: $M_b^{-1}=\tfrac{1}{1+m_P}(D_{ov}^\dagger+\nu)^{-1}$, so $[\![M_b^{-1}]\!]=\tfrac{1}{2(1+m_P)}$; and
  $[\![(1-D_{ov}^\dagger)M_b^{-1}]\!]=\tfrac{1}{1+m_P}\cdot\tfrac{\nu-1}{2}=-\tfrac{1+2m_P}{2(1+m_P)^2}$
  (using $\nu-1=-\tfrac{1+2m_P}{1+m_P}$).

### Channel contacts ($/V_{st}$) and csub $=-$(contact)
`etadag_xi`$=-1$ (both runs). `xidag_eta`$=-2[\![M_b^{-1}]\!]$. `xidag_1mDdag_eta`$=-2[\![(1-D^\dagger)M_b^{-1}]\!]$.
$$
\text{csub}_{PS}=1+2[\![M_b^{-1}]\!],\qquad \text{csub}_{FS}=1-2[\![(1-D_{ov}^\dagger)M_b^{-1}]\!].
$$

| | csub$_{PS}$ | csub$_{FS}$ |
|---|---|---|
| mF | $2$ | $2-m_F$ |
| mP | $1+\tfrac{1}{1+m_P}=\tfrac{2+m_P}{1+m_P}$ | $1+\tfrac{1+2m_P}{(1+m_P)^2}=2-\left(\tfrac{m_P}{1+m_P}\right)^2$ |

### Numerics ($m_F=0.1$, $m_P=0.1i$)
- csub$_{PS}$(mF) $=2$;  csub$_{FS}$(mF) $=1.9$.
- csub$_{PS}$(mP) $=\tfrac{2+0.1i}{1+0.1i}=1.990099-0.099010\,i$.
- csub$_{FS}$(mP) $=2-\left(\tfrac{0.1i}{1+0.1i}\right)^2=2-(-0.009705+0.0019606\,i)=2.009705-0.0019606\,i$.

### Notebook fix (`jj_validate_{mF,mP}` cond cells)
Correct already: csub$_{PS}$(mF)$=2$, csub$_{PS}$(mP)$=1+\tfrac{1}{1+m_P}$. WRONG: csub$_{FS}$(mF) -> use
$2-m_F$ (not $1+\tfrac{1}{1+m_F}$); csub$_{FS}$(mP) -> use $2-(\tfrac{m_P}{1+m_P})^2$ (not $2$). After the fix,
each residual is the pure dynamical part (e.g. mF-FS $\to m_F\delta T=-0.0038$).

## 11. Why mF-FS is NOT zero -- the parity Ward identity (anomaly), and how to get the true zero
The intuition "$\langle\sigma_{FS}\rangle=0$ in the parity-invariant mF run" is right IN THE CONTINUUM, where
$\sigma_{FS}\to\eta^\dagger\xi-\xi^\dagger\eta$ is PURELY parity-odd ($\to-\sigma_{FS}$ under parity). On the
LATTICE $\sigma_{FS}$ is NOT purely parity-odd, because of the $(1-D_{ov})$ dressing.

### Parity transform of $\sigma_{FS}$ (exact)
Parity (eq:parity_Psi/Psibar): $U\equiv\sigma_1 P$, $U^2=1$; $\xi\to U\eta,\ \eta\to U\xi,\ \xi^\dagger\to\eta^\dagger U,\ \eta^\dagger\to\xi^\dagger U$;
operator condition $U D_{ov}U=D_{ov}^\dagger$ (eq:parity_cond_D), so $U D_{ov}^\dagger U=D_{ov}$. Then
$$
\eta^\dagger\xi\to\xi^\dagger\eta,\qquad
\xi^\dagger(1-D_{ov}^\dagger)\eta\to\eta^\dagger U(1-D_{ov}^\dagger)U\xi=\eta^\dagger(1-D_{ov})\xi,
$$
$$
\Rightarrow\ \sigma_{FS}^{P}=\xi^\dagger\eta-\eta^\dagger(1-D_{ov})\xi
   = -\sigma_{FS}+\underbrace{(\eta^\dagger D_{ov}\xi+\xi^\dagger D_{ov}^\dagger\eta)}_{=\,S_2\ (\text{massless action})}.
$$
So $\boxed{\sigma_{FS}^{P}=-\sigma_{FS}+S_2}$: parity-odd PLUS the massless action. The "$+S_2$" is the GW/lattice
parity-breaking term ($\to0$ in the continuum).

### Ward identity
The mF action $S_2+m_F\sigma_{PS}$ is parity-invariant, so $\langle\sigma_{FS}\rangle=\langle\sigma_{FS}^P\rangle=-\langle\sigma_{FS}\rangle+\langle S_2\rangle$:
$$
\boxed{\ \langle\sigma_{FS}\rangle_{mF}=\tfrac12\langle S_2\rangle\ }\quad(\text{exact, to Zolotarev}).
$$
Not zero. Direct evaluation: $\langle S_2\rangle=-2N+2m_F T_R$ (with $T_R=\operatorname{Re}\operatorname{tr}[A D_{m_F}^{-1}]$,
$N=2V_{st}$), so $\langle\sigma_{FS}\rangle=m_F T_R-N$ -- matches §8. Using $\sigma_{PS}=-2T_R$ this is EXACTLY the
§8 relation $\langle\sigma_{FS}\rangle=-\tfrac{m_F}{2}\sigma_{PS}-2V_{st}$. So that relation IS the parity Ward
identity, and the post-csub residual $m_F\delta T$ is the non-contact part of the parity-even anomaly
$\tfrac12\langle S_2\rangle$ -- a finite-$a$ lattice artifact, NOT Zolotarev, NOT the order parameter.

### The genuinely zero (parity-odd) order parameter
Parity-odd projection $O_{\rm odd}=\tfrac12(O-O^P)$: $\sigma_{FS}-\tfrac12 S_2=\sigma_{FS}+\tfrac{m_F}{2}\sigma_{PS}+2V_{st}$.
By parity its VEV is EXACTLY zero. Numerically (dense, $m_F=0.1$):
$-1.9038+0.05(-1.9238)+2=+1\times10^{-5}$ (Zolotarev floor). So PLOT
$$
\langle\sigma_{FS}\rangle+\tfrac{m_F}{2}\langle\sigma_{PS}\rangle+2V_{st}\ \xrightarrow{}\ 0
$$
to confirm "no parity-odd condensate" -- the symmetry statement you expected. The contact-only csub leaves
$m_F\delta T$ because it does NOT remove the dynamical half of $\tfrac12\langle S_2\rangle$ ($=-\tfrac{m_F}{2}\sigma_{PS}$);
subtract the FULL $\tfrac12 S_2$ to recover the protected zero.

### Why the residual is $O(m^2)$
Write the Ward identity via $\langle S_2\rangle=-2N-m_F\langle\sigma_{PS}\rangle$ ($\langle S_{tot}\rangle=-2N=-$(total
$(\xi,\eta)$ dof), $S_2=S_{tot}-m_F\sigma_{PS}$; $N=\operatorname{tr}[A\mathbb 1]=2V_{st}$):
$$
\langle\sigma_{FS}\rangle=-N-\tfrac{m_F}{2}\langle\sigma_{PS}\rangle,\qquad
\langle\sigma_{PS}\rangle=\underbrace{-N}_{\text{contact}\,(-2V_{st})}+\underbrace{(\sigma_{PS}\text{ shift})}_{=-2\delta T}.
$$
After removing the (mass-dependent) contact, the residual FACTORIZES:
$$
\sigma_{FS}^{\rm sub}=-\tfrac{m_F}{2}\times(\sigma_{PS}\text{ shift})=O(m^2),
$$
two powers of $m$: (1) $-\tfrac{m_F}{2}$ = the parity-ANOMALY coupling (the cross channel sees the mass only
through the $-m_F\langle\sigma_{PS}\rangle$ term in $\langle S_2\rangle$); (2) the $\sigma_{PS}$ shift $=-2\delta T=O(m_F)$
= the conjugate channel's LINEAR response to its own source ($\delta T=-m_F\operatorname{Re}\operatorname{tr}[A D_{ov}^{-2}]+\dots$,
vanishes at $m=0$). General principle: a WRONG-symmetry-channel response is one power of $m$ below the direct
(susceptibility) response. Same statement as the mP real part: parity-even (wrong-symmetry) -> no $O(m_P)$
term (that is the odd order parameter) -> starts at $O(m_P^2)$.

## 12. What changes WITH SSB (the whole point)
Everything in §1-11 is the SYMMETRIC-phase (free) baseline: UV contacts $+$ $O(m)$ susceptibility responses
that vanish as $m\to0$. SSB is an IR effect; the UV/contact structure (and hence the csub) is UNCHANGED
($\tfrac12$ contact is GW, gauge-/phase-independent, §6a). What changes:

1. **Order-parameter PLATEAU (the signature).** In the SYMMETRIC channel the response stops vanishing:
   $\langle\sigma\rangle^{\rm sub}(m)\to\Sigma\neq0$ as $m\to0^+$ (Banks-Casher), instead of $O(m)\to0$. A
   discontinuity at $m=0$ ($\langle\sigma\rangle^{\rm sub}(0^+)-(0^-)=2\Sigma$); the susceptibility
   $d\langle\sigma\rangle/dm$ DIVERGES (the linear slope becomes a step). Order of limits: $V\to\infty$ FIRST,
   then $m\to0$ (finite $V$ has no true SSB -> need finite-size scaling, a la Karthik-Narayanan).
2. **Spectral picture (§6c).** $\Sigma=\pi\rho(0)$: SSB <=> near-zero modes of $D_{ov}$ accumulate,
   $\rho(0)\neq0$. The free theory is gapped ($\rho(0)=0$) -> no SSB -> the plateau is absent (what the free
   tests confirm). $\delta T$ stops being $\propto m$ and tends to a constant.
3. **Cross-channel power FLIPS $m^2\to m$.** The parity Ward identity (§11) is EXACT with interactions, so
   cross residual $=-\tfrac{m}{2}\langle\sigma_{\rm conj}\rangle^{\rm sub}$. Symmetric: $\langle\sigma_{\rm conj}\rangle^{\rm sub}=O(m)$
   -> cross $=O(m^2)$. Broken: $\to\Sigma$ -> cross $=-\tfrac{m}{2}\Sigma=O(m)$. So the cross-channel SCALING is
   itself an SSB diagnostic (it still $\to0$ as $m\to0$ -- the wrong-symmetry order parameter stays zero --
   but linearly, not quadratically). The protected zero $\sigma_{FS}-\tfrac12 S_2=0$ still holds (parity
   unbroken by a parity-even mass).
4. **Goldstones (continuous breaking).** Flavor $SU(2)\to U(1)$ by $\langle\sigma_{PS}\rangle$ -> Goldstone
   modes; disconnected susceptibility $\chi_{\rm disc}\sim\Sigma^2/m$ DIVERGES as $m\to0$; large IR
   fluctuations of the condensate (gauge noise, not stochastic -- §condensate_impl_plan caveat).
5. **QED3 specifics.** mF run probes FLAVOR SSB (does $\langle\sigma_{PS}\rangle$ condense?) -- Karthik-Narayanan
   (arXiv:1512.02993) find NO bilinear condensate for small $N_f$ ($N_f=2$ conformal), so expect the plateau
   to be ABSENT and $\langle\sigma_{PS}\rangle^{\rm sub}\to0$. mP run probes PARITY: parity is EXPLICITLY broken
   at finite $a$ by the overlap anomaly (the $S_2$/Chern-Simons term, §11), so genuine parity SSB must be
   disentangled from that lattice/anomaly piece via $a\to0$ -- the imaginary $\langle\sigma_{FS}\rangle$
   extrapolated $m_P\to0$ then $a\to0$.
