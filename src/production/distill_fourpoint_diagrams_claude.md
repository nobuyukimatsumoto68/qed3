# The ten $\sigma\sigma$ four-point diagrams as distillation traces (chunk 4)

Translation of the ten diagrams A--J of `qed3int_v3-4.pdf` Ch.5 (Eqs. 5.4--5.20, Fig. 1) into exact
distillation traces over the forward perambulator $\tau$, the furnished perambulator $\tau'$, and the
$\ell{=}0$ elemental $\Phi$. This is the plan for chunk 4 (`distill_contract_claude.py` four-point +
coupled GEVP). It supersedes Section 5 of `distillation_for_two_meson_claude.md`, which used a BACKWARD
perambulator $\bar\tau$ -- see the correction below.

## 0. Objects (from chunks 1-3, all $N_v\times N_v$ per config, spin folded into the mode index)
- $\tau(t',t)=V^\dagger(t')D_{ov}^{-1}V(t)$ -- forward perambulator (stored, windowed).
- $\tau'(t',t)=V^\dagger(t')(1-D_{ov}^\dagger)D_{ov}^{-1}V(t)$ -- furnished ("gw") perambulator (stored).
- $\Phi(t)=V^\dagger(t)\,W_0\,V(t)$, $W_0=\mathrm{diag}_x(A_x)\otimes\mathbb 1_{\rm spin}$, $A_x=\text{dual\_area}_x\,Y_{00}$
  -- the $\ell{=}0$ area-weighted vertex (one per vertex; higher $\ell$ carry $Y_{\ell m}$).

## 0b. What $\tau'$ is (and is NOT) -- the furnished perambulator

$\tau'$ is **not** a time-reversed perambulator. Time reversal of $\tau(t',t)$ would swap source/sink,
$\tau(t,t')$, or reflect $t\to N_t-t$; $\tau'$ is a **different object at the SAME $(t',t)$**. It is the
**FS-furnished** perambulator: the FS scalar vertex is $\tilde S_{\rm FS}=-(1-D_{ov}^\dagger)$ (Eq. 5.3), so
wherever an FS vertex sits it dresses the adjacent (sink-side) leg with the operator $(1-D_{ov}^\dagger)$.
Side by side:
$$
\tau(t',t) = V^\dagger(t')\,D_{ov}^{-1}\,V(t),
\qquad
\tau'(t',t) = V^\dagger(t')\,(1-D_{ov}^\dagger)\,D_{ov}^{-1}\,V(t).
$$

**Computation (cheap, no new solve).** Both come from the SAME forward solve $\psi_l=D_{ov}^{-1}V(t)$; then
$$
\tau(t',t)=V^\dagger(t')\,\psi\big|_{t'},
\qquad
\tau'(t',t)=V^\dagger(t')\,\big[(1-D_{ov}^\dagger)\psi\big]_{t'} ,
$$
i.e. $\tau'$ is $\tau$ with ONE extra application of $(1-D_{ov}^\dagger)$ at the sink -- exactly the
`tau_gw` computed and stored in chunk 2 (`distill_peram_claude.cu`).

**Formal relation to $\tau$.** Using the GW/overlap identity (IV.17) $D_{ov}^{-\dagger}=1-D_{ov}^{-1}$, multiply
by $D_{ov}^\dagger$ on the left: $1=D_{ov}^\dagger-D_{ov}^\dagger D_{ov}^{-1}$, hence $D_{ov}^\dagger D_{ov}^{-1}=D_{ov}^\dagger-1$, so
$$
(1-D_{ov}^\dagger)\,D_{ov}^{-1} = D_{ov}^{-1}-D_{ov}^\dagger D_{ov}^{-1} = 1 + D_{ov}^{-1} - D_{ov}^\dagger ,
$$
$$
\boxed{\ \tau'(t',t) = \tau(t',t) + O(t',t) - \Xi(t',t)\ },\qquad
O(t',t)\equiv V^\dagger(t')V(t),\quad \Xi(t',t)\equiv V^\dagger(t')D_{ov}^\dagger V(t).
$$
$O$ is the mode-overlap ($O(t,t)=\mathbb 1$, since $V^\dagger(t)V(t)=\mathbb 1$; off-diagonal in $t$ it is not
$\mathbb 1$). $\Xi$ is one $D_{ov}^\dagger$ matrix element -- but $D_{ov}^\dagger$ is NONLOCAL (the Zolotarev
sign function), so $\tau'$ does **not** reduce to $\tau$; it is a genuine separate object (hence we store both).

**Note on the (IV.17) "collapse".** The clean collapse $(1-D_{ov}^\dagger)D_{ov}^{-\dagger}=-D_{ov}^{-1}$ (plan
key-simpl #2) is for the BACKWARD inverse $D_{ov}^{-\dagger}$. The four-point uses only the FORWARD
$D_{ov}^{-1}$ (Sec. 1a), so that collapse never triggers here -- FS legs always use $\tau'$ directly.

**No $\Box$ on the $\psi$-side for the furnishing.** $(1-D_{ov}^\dagger)$ acts on the EXACT solved propagator
$\psi=D_{ov}^{-1}V(t)$, projected once by $V^\dagger(t')$ -- there is NO intermediate $\Box=VV^\dagger$ between
$D_{ov}^{-1}$ and $(1-D_{ov}^\dagger)$. This is forced: $(1-D_{ov}^\dagger)$ is time-nonlocal, so it cannot be a
single-timeslice elemental $V^\dagger(t)(1-D_{ov}^\dagger)V(t)$ -- it must dress the propagator (hence $\tau'$).
The LOCAL area $W_0$ is different: it stays a smeared elemental $\Phi=V^\dagger W_0V$ (its $\Box$ is the
intended operator smearing). At an FS vertex the product $\Phi(x)\tau'(x,\cdot)$ thus reinserts $\Box(x)$
between $W_0$ and $(1-D_{ov}^\dagger)$ -- EXACT at $N_v{=}2N_s$ (L1), the standard smearing approximation at
truncated $N_v$. (Optional L2+ refinement: fold $W_0(1-D_{ov}^\dagger)$ into ONE sink operator on $\psi$ to
drop that $\Box$ too; identical at L1.)

**Guiding principle (generalized distillation for nonlocal $\Gamma$).** This is NOT vanilla Peardon-2009
distillation (which assumes LOCAL $\Gamma$ -> elemental $\Phi=V^\dagger\Gamma V$). A nonlocal $\Gamma$ (here
$1-D_{ov}^\dagger$) is instead folded into a GENERALIZED perambulator (dress the exact propagator, smear
$\Box$ only at the ends). Used CONSISTENTLY -- same prescription for every nonlocal $\Gamma$ and every local
vertex, in all diagonal AND off-diagonal GEVP correlators -- it is a valid variational basis: GEVP
eigenvalues (masses) and the light branch's $F^2$/$\sigma\sigma$ content are physical (smearing changes
overlaps, not states), and at $N_v{=}2N_s$ (L1) it is exact.

## 1. TWO corrections vs the old method note

**(a) Forward-only strategy -- NO backward perambulator.**

*The problem.* A meson/four-point loop runs in BOTH time directions. The single meson $C_S$ is a loop
$0\to t\to 0$:
$$
C_S(0{\to}t)=\mathrm{Tr}\big[\Phi(0)\,\underbrace{\tau(0,t)}_{t\to 0\ \text{leg}}\,\Phi(t)\,\underbrace{\tau(t,0)}_{0\to t\ \text{leg}}\big].
$$
In a normal (Wilson / point-source) calculation one computes $G(\cdot;0)$ from a source at $0$ and gets the
RETURN leg $G(0;t)$ for free via $\gamma_5$-hermiticity, $G(0;t)=\gamma_5 G(t;0)^\dagger\gamma_5$. **We have
no $\gamma_5$-hermiticity** (2-component GW), so naively the return leg would need a BACKWARD inverse
$D_{ov}^{-\dagger}$ (the $\bar\tau$ of the old method note).

*The distillation resolution.* The perambulator is built by SOLVING AT EVERY SOURCE TIMESLICE in the window,
$\tau(t',t)=V^\dagger(t')D_{ov}^{-1}V(t)$ for all $t,t'\in[t_0,t_0{+}W)$. So both directions are already
stored, BOTH as FORWARD inverses:
- $\tau(t,0)$ = forward $D_{ov}^{-1}$, source at $0$, sink read at $t$;
- $\tau(0,t)$ = forward $D_{ov}^{-1}$, source at **$t$**, sink read at $0$.
The return leg $\tau(0,t)$ is NOT backward -- it is a plain forward solve whose source sits on the later
timeslice. Because we solved at all window sources, we have it. **$D_{ov}^{-\dagger}$ appears nowhere.**

*Why every diagram is forward.* Ch.5 writes all ten diagrams (5.7--5.20) with $D_{ov}^{-1}$ only (traces
$\mathrm{tr}\,SG\cdots SG$, $G=D_{ov}^{-1}$); each timeslice projector $\Pi_s$ just says "this leg's
source/sink is $s$", i.e. selects the stored forward $\tau(\text{sink},\text{source})$. So
$S_S,T_S,V_S,C_S,\dots$ are products of $\Phi$'s and forward $\tau$'s -- exactly the code
(`M = Phi0 @ tau[0,dt] @ Phit @ tau[dt,0]` for $C_S$'s two legs).

*$\tau'$ and (IV.17).* The FS furnished $\tau'=V^\dagger(1-D_{ov}^\dagger)D_{ov}^{-1}V$ is ALSO forward (a
forward inverse with the vertex operator at the sink). The backward collapse
$(1-D_{ov}^\dagger)D_{ov}^{-\dagger}=-D_{ov}^{-1}$ (IV.17) would matter only if a backward leg appeared -- it
never does -- so (IV.17) is used only as the T2b validation check, NOT in the contraction.

*The price.* Forward-only costs the all-to-all: $W\times N_v$ solves per config (solve at every window
source) instead of one. That is the trade for having no $\gamma_5$-herm, and it is cheap here -- the L1
overlap solves are small and the perambulator is reused across ALL diagrams and the whole GEVP.

**(b) FS furnishing = the $\tau'$ leg, one per $\tilde S$ vertex.** $G_4=\langle S_4\rangle+\langle\tilde S_4\rangle$
(5.5): $S_4$ from the $\eta^\dagger S\xi$ half of each $\sigma$ (vertex $S=1$), $\tilde S_4$ from the
$\xi^\dagger\tilde S\eta$ half (vertex $\tilde S$). For a leg entering a vertex $v$:
- $S=1$ (both PS and FS $S_4$; and PS $\tilde S_4$): plain leg $\tau(v,\cdot)$, vertex $\Phi(v)$.
- $\tilde S=-(1-D_{ov}^\dagger)$ (FS $\tilde S_4$): the furnishing sits at the vertex $v$ (sink), so the
  incoming leg is the furnished $-\tau'(v,\cdot)$, and the vertex still carries $\Phi(v)$ (the area $A_v$).
  [At $N_v=2N_s$ the $VV^\dagger=\mathbb 1$ insertion that separates $\Phi$ from $\tau'$ is EXACT.]

So: **$S_4$** = the ten diagrams with all legs $\tau$; **$\tilde S_4^{\rm PS}$** = identical (all $\tau$);
**$\tilde S_4^{\rm FS}$** = the ten diagrams with all legs $-\tau'$ (each vertex furnished). Hence
$$
G_4^{\rm PS\cdot PS}=2\,G_{10}[\tau],\qquad
G_4^{\rm FS\cdot FS}=G_{10}[\tau]+G_{10}[-\tau'],
$$
where $G_{10}[\,\cdot\,]$ is the weighted ten-diagram sum with the indicated leg object. (Mixed FS$\cdot$PS
and the $F^2$ coupling: Section 4.)

## 2. Vertex / area rule (from `sigma_sigma_connected_vertex_note_claude.md`)
\#($A$ factors) = \#($\sum_x$) = \#(vertices) = one $\Phi$ per vertex. A connected loop through $n$ vertices
has $n$ $\Phi$'s and $n$ legs; a disconnected product carries its $\Phi$'s split across independent factors.
Distillation is EXACT, so squares/products use distinct exact factors (no distinct-hit trick, no square bias).

## 3. The ten diagrams (weights $w_i$ from Fig. 1)

Building blocks (5.15--5.20), plain-$\tau$ (i.e. $S$ / PS) versions:
$$
D_S(t)=\mathrm{Tr}[\Phi(t)\tau(t,t)],\qquad
D'_S(t)=\mathrm{Tr}[\Phi(t)\tau(t,t)\Phi(t)\tau(t,t)],
$$
$$
C_S(0{\to}t)=\mathrm{Tr}[\Phi(0)\tau(0,t)\Phi(t)\tau(t,0)],\qquad
V_S(0{\to}t)=\mathrm{Tr}[\Phi(0)\tau(0,0)\Phi(0)\tau(0,t)\Phi(t)\tau(t,0)],
$$
$$
S_S(0{\to}t)=\mathrm{Tr}[\Phi(0)\tau(0,0)\Phi(0)\tau(0,t)\Phi(t)\tau(t,t)\Phi(t)\tau(t,0)],
$$
$$
T_S(0{\to}t)=\mathrm{Tr}\big[\big(\Phi(0)\tau(0,t)\Phi(t)\tau(t,0)\big)^2\big].
$$

| # | $w$ | Ch.5 | distillation (plain $\tau$) | type |
|---|---|---|---|---|
| A | 4 | 5.7  $=-S_S$ | $-S_S(0{\to}t)$ | connected, nested |
| B | 2 | 5.8  $=-T_S$ | $-T_S(0{\to}t)$ | connected, nested |
| C | 4 | 5.9  | $D_S(t)\,V_S(0{\to}t)$ | loop $\times$ nested |
| D | 4 | 5.9' | $D_S(0)\,V_S(t{\to}0)$ | time-reverse of C |
| E | 2 | 5.10 | $C_S(0{\to}t)^2$ | factorizable |
| F | 1 | 5.11 | $D'_S(0)\,D'_S(t)$ | factorizable |
| G | 4 | 5.12 | $-\,D_S(0)\,D_S(t)\,C_S(0{\to}t)$ | factorizable |
| H | 1 | 5.13 | $-\,D_S(0)^2\,D'_S(t)$ | factorizable |
| I | 1 | 5.13' | $-\,D_S(t)^2\,D'_S(0)$ | time-reverse of H |
| J | 1 | 5.14 | $D_S(0)^2\,D_S(t)^2$ | factorizable |

$S_4(0,t)=\sum_i w_i\mathfrak A_i$. The FS $\tilde S_4$ uses the SAME table with every $\tau\to-\tau'$ (so
each building block's legs become furnished: e.g. $C_S^{\rm FS}=\mathrm{Tr}[\Phi(0)\tau'(0,t)\Phi(t)\tau'(t,0)]$,
$D'^{\rm FS}_S(t)=\mathrm{Tr}[\Phi\tau'(t,t)\Phi\tau'(t,t)]$, sign $(-1)^{\#\rm legs}$ per term).

## 4. Coupled GEVP basis $\{F^2,\ \sigma_{\rm FS}\sigma_{\rm FS},\ \sigma_{\rm PS}\sigma_{\rm PS},\ \sigma_{\rm FS}\sigma_{\rm PS}\}$
- $\langle\sigma_a\sigma_a(t)\,\sigma_b\sigma_b(0)\rangle$: for each (a,b) channel, the ten-diagram sum with the
  appropriate leg objects per timeslice-half (PS-half $\to\tau$, FS-half $\to-\tau'$). PS$\cdot$PS $=2G_{10}[\tau]$;
  FS$\cdot$FS $=G_{10}[\tau]+G_{10}[-\tau']$; the mixed blocks pick $\tau$ on the PS operator's legs and
  $-\tau'$ on the FS operator's legs (to be written explicitly per diagram in code).
- $\langle F^2(t)\,\sigma_a\sigma_a(0)\rangle$: the MIXING cross -- gluonic $O_F(t)$ times the equal-time
  scalar $[D_S(0)^2+D'_S(0)]$ (PS: $\tau$; FS: $\tau'$), exactly the chunk-3 loops (now exact, de-noised).
- $\langle F^2(t)F^2(0)\rangle$: the standalone glue $F^2$ correlator (existing `glue_f2_v2_shapes`).
Solve the generalized eigenvalue problem $C(t)v=\lambda(t)C(t_0)v$; the light branch's overlap onto both
$F^2$ and $\sigma\sigma$ is the mixing; $\Delta_-$ from its plateau.

## 5. Validation (chunk 4)
- **T4a (factorizable cross-checks):** the distilled $E,F,G,H,I,J$ built from the FULL contraction must equal
  their loop-built products ($E=C_S^2$, $F=D'_S D'_S$, $J=D_S^2 D_S^2$, ...) to machine precision -- pins the
  index/spin bookkeeping. (Trivially true by construction here, so the real test is that the SAME $\Phi,\tau$
  reproduce the chunk-3 loops AND assemble into A/B.)
- **T4b:** time-reversal $D(t)=C(N_t{-}t)$, $I(t)=H(N_t{-}t)$ numerically.
- **T4c:** GEVP matrix Hermitian/PSD; its $F^2$ diagonal block equals the standalone glue $F^2$.
- **X-check:** free field (U=1) vs the analytic free four-point where available.

## 6. TO CONFIRM with NM before coding
1. **A$=-S_S$, B$=-T_S$** traces (Sec. 3) -- sign and $\Pi$ ordering read from (5.7)/(5.8),(5.19)/(5.20).
2. **C,D** $=D_S\cdot V_S$ (a disconnected loop $\times$ a nested 3-leg loop, gauge-connected) -- confirm this
   is the intended factorization of (5.9), not a single connected 4-leg object.
3. **FS furnishing placement:** $\tilde S=-(1-D_{ov}^\dagger)$ enters as the incoming leg $-\tau'(v,\cdot)$ at
   each $\tilde S$ vertex, with $\Phi(v)$ still carrying $A_v$ (Sec. 1b). Confirm vs the FS single-meson
   $C_S^{\rm FS}=\mathrm{Tr}[\Phi\tau'\Phi\tau']$ (both vertices furnished, as in the existing $V^{--}$ conn).
4. **Signs** on the FS terms $(-1)^{\#\rm furnished\ legs}$ and the overall diagram signs (A,B,G,H,I negative).
