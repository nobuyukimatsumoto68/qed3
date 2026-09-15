# FS-channel $\sigma^2$ two-point, diagram-by-diagram -- what changes vs PS

> **!!! SUPERSEDED / WRONG (2026-09-09b). This note's derivation below uses the FORWARD furnished leg
> $-\tau'=-(1-D_{ov}^\dagger)D^{-1}$ (`taugw`) -- that is INCORRECT.** It gives a spurious FS one-point
> $0.549$ (should be $0$), a live tadpole, and a huge diagram G. The correct FS uses the **BACKWARD inverse**
> $D^{-\dagger}=1-D^{-1}=O-\tau$; by GW $(1-D_{ov}^\dagger)D^{-\dagger}=-D^{-1}=-\tau$ EXACTLY (no new solve).
> The correct rule: **the $S/\tilde S$ kernels sum PER CLOSED LOOP** (not per diagram) and **always use the
> $\tfrac12$-improved propagator** $\tilde\tau=\tau-\tfrac12 I$; then each loop $\to 2\times$ and every TADPOLE
> loop vanishes, so diagrams C,D,G,H,I,J $\to 0$ and only A,B,E(+F vacuum) survive. Implementation:
> FS diagram$_i = 2^{n_\text{loops}[i]}\times$`diags_pair(tau, CONTACT=0.5)`. See `fs_diag_corr_v2_claude.py`
> (correct diagram-by-diagram + triple subtraction), `fs_channels_v2_claude.py` (corrected 3-op channels/GEVP),
> and the memory `project_sigma_sigma_f2_mixing.md` "SESSION 2026-09-09 (b)". The material below is kept only
> for the operator definitions and the (still-valid) $S/\tilde S$ non-mixing bookkeeping; ignore its forward
> `taugw` legs and its "furnished $\tfrac12$" contact.

Following `sigma2_diagram_vacsub_method_claude.md` (step 1 = per-diagram linear, plateau-subtracted), but
the FS channel differs from PS in two ways. One is settled (I will implement it); one is a genuine open
choice for you.

## Difference 1 (SETTLED by the doc -- I will implement): FS$\cdot$FS is a SUM of two $G_{10}$'s

$\sigma=\eta^\dagger S\xi+\xi^\dagger\tilde S\eta$, with $S=1$ and $\tilde S=1$ (PS) / $\tilde S=-(1-D_{ov}^\dagger)$
(FS). The four-point splits as $G_4=\langle S_4\rangle+\langle\tilde S_4\rangle$ (Ch.5 Eq 5.5) with **no cross
terms** ($S$ and $\tilde S$ never mix in a trace). So per diagram $i$ (weights $W_{10}$ unchanged):
$$
\text{PS}\cdot\text{PS} = 2\sum_i W_{10}[i]\,G_{10}[\tau]_i,\qquad
\text{FS}\cdot\text{FS} = \sum_i W_{10}[i]\big(G_{10}[\tau]_i + G_{10}[-\tau']_i\big),
$$
where $\tau=$ forward perambulator (`tau`), $\tau'=$ furnished perambulator (`tau_gw`), and $G_{10}[\text{leg}]_i
= $ `diags_pair(Phi, leg, s, t)[i]`. **The per-diagram FS correlator is $G_{10}[\tau]_i+G_{10}[-\tau']_i$**
(the PS diagram plus its $-\tau'$ partner), NOT $2\,G_{10}[-\tau']_i$.

> **BUG in the existing tools for FS.** `diag_effmass_claude.py` / `diag_corr_linear_claude.py` set
> `leg=-taugw` and return `2 * sum W10 * diags_pair(-taugw)` $=2\,G_{10}[-\tau']$ -- this DROPS the
> $G_{10}[\tau]$ half and double-counts $-\tau'$. The FS driver must sum the two legs. (PS is correct:
> $S=\tilde S=1$ makes both halves $G_{10}[\tau]$, giving the factor 2.)

## Difference 2 (OPEN -- your call): the equal-time CONTACT for the $\tau'$ leg

PS subtracts the GW contact $\tfrac12$ on every equal-time diagonal (`tt=tau[a,a]-0.5 I`), which makes the
PS tadpole $\langle D_S^{PS}\rangle=0$ exactly (contact-saturated). For the **furnished $\tau'$ leg the
equal-time contact is NOT $\tfrac12$** (memory: "FS contact unsolved; forward $\tau'$ contact $\ne1/2$ and is
non-scalar"). The $S$-part ($\tau$ leg) keeps $\tfrac12$; the question is what to subtract on the $-\tau'$
diagonal $-\tau'(s,s)$:

**DECISION (NM 2026-09-09): subtract the overlap $\tfrac12$ contact on the $\tau'$ leg too -- BUT the
furnished factor $(1-D_{ov}^\dagger)$ multiplies it, so the FS subtraction is a full matrix, not $\tfrac12 I$.**

Derivation. The furnished perambulator is
$$
\tau' = V^\dagger\,(1-D_{ov}^\dagger)\,D_{ov}^{-1}\,V \quad(\texttt{tau\_gw}).
$$
The overlap propagator has the ultralocal GW contact $D_{ov}^{-1}(x,y)=\tfrac12\delta_{xy}+G_{\rm reg}(x,y)$
(from $D_{ov}^{-1}+D_{ov}^{-\dagger}=1$; $\mathrm{Re}\,\mathrm{diag}\,D_{ov}^{-1}=\tfrac12$ EXACT). Because the
$\tfrac12\delta$ sits to the RIGHT of the furnished factor, the coincident (contact) piece of $\tau'$ is that
factor times $\tfrac12$:
$$
\mathrm{contact}(\tau') = \tfrac12\,V^\dagger(1-D_{ov}^\dagger)V = \tfrac12\,(O-\Xi),
\qquad O=V^\dagger V,\ \ \Xi=V^\dagger D_{ov}^\dagger V .
$$
Using the GW identity $D_{ov}^\dagger D_{ov}^{-1}=D_{ov}^\dagger-1$ one gets the exact (any-$L$) relation
$\tau' = \tau + O - \Xi$, i.e. $O-\Xi=\tau'-\tau$. Hence
$$
\boxed{\ \mathrm{contact}(\tau') = \tfrac12\,(\tau'-\tau)\ }\qquad\text{(computable from the stored \texttt{tau}, \texttt{tau\_gw}).}
$$
The FS $\tilde S$-part leg is $-\tau'$, so its contact is $-\tfrac12(\tau'-\tau)$; normal-ordering (removing it)
gives the equal-time diagonal
$$
\boxed{\ \tilde\tau^{\,\tilde S}(s,s) = -\tau'(s,s) + \tfrac12\big(\tau'(s,s)-\tau(s,s)\big)
      = -\tfrac12\big(\tau'(s,s)+\tau(s,s)\big) = -\tfrac12\big(\texttt{tau\_gw}[s,s]+\texttt{tau}[s,s]\big).\ }
$$
Off-diagonal legs are untouched: $-\tau'(s,t)=-\texttt{tau\_gw}[s,t]$. The $S$-part leg is $\tau$ with the
ordinary contact $\tilde\tau^{S}(s,s)=\tau(s,s)-\tfrac12 O$ (at L1 $O=V^\dagger V=I$, so $\tau(s,s)-\tfrac12 I$;
the existing PS code is the L1 form).

So the FS driver does NOT call `diags_pair(-taugw, CONTACT=0.5)` (that would wrongly use $-\tau'(s,s)-\tfrac12 I$).
It uses a custom $\tilde S$ diagram builder with equal-time diagonal $-\tfrac12(\texttt{tau\_gw}+\texttt{tau})[s,s]$
and off-diagonals $-\texttt{tau\_gw}[s,t]$; the $S$-part is the ordinary `diags_pair(tau, CONTACT=0.5)`.

DIAGNOSTIC: print the FS tadpole $\langle\mathrm{Tr}[\Phi\,\tilde\tau]\rangle$ for each part -- for the $S$-part
($\tau$) it is $0$; for the $\tilde S$-part it should be small if the furnished contact is correctly removed.

The vacuum-free diagram SELECTION (keep A,B,C,D,E,G; drop F,H,I,J) is unchanged from PS -- the contact only
sets where the tadpole-carrying diagrams sit.

## Note on the $\tilde S=-(1-D_{ov}^\dagger)$ overall MINUS sign (innocuous for the two-point)

$\tilde S=-(1-D_{ov}^\dagger)$ puts one factor $(-1)$ on each $\sigma_{FS}$ vertex. **Every diagram of
$\langle\sigma^2\sigma^2\rangle$ has exactly 4 $\sigma$ vertices (4 propagator legs)**, so the $\tilde S$-part
carries $(-1)^4=+1$: the overall sign cancels. VERIFIED numerically (2026-09-09): computing each of the 10
diagrams with leg $=-\tau'$ vs $+\tau'$ (matched contact $\mp\tfrac12(\tau'+\tau)$, $\mp\tau'(s,t)$) gives
IDENTICAL results (diff $=0$ to machine precision) for A..J. So `leg=-taugw` is correct and the result is
sign-independent for this two-point. The minus would matter only for an ODD-vertex correlator; the FS tadpole
$D_S^{FS}$ (one vertex, sign-dependent) enters the diagrams only in even combinations ($D_S^2$ in G,H,I,J;
$D_S\!\cdot\!V_S$ with 4 total legs in C,D), so the physics is unaffected.

## Cross $\langle O_{\sigma\sigma}^{FS}\,\sigma_{FS}^2\rangle$ (FS two-$\sigma$; NM 2026-09-09)

$O_{\sigma\sigma}^{FS}(t)=\sigma_{FS,00}(t)\,\sigma_{FS,00}(t{+}\delta)$ (time-split, $\delta=$ `SPLIT`, default 1) --
the genuine four-fermion two-meson interpolator in the FS channel. The cross with $\sigma_{FS}^2(s)$ has 4
$\sigma_{FS}$ vertices, so by the $S/\tilde S$-non-mixing rule it splits into $\langle S_4\rangle+\langle\tilde S_4\rangle$
= the all-$\tau$ two-meson plus the all-$(-\tau')$ two-meson:
$$
\langle O_{\sigma\sigma}^{FS}(t)\,\sigma_{FS}^2(s)\rangle_c
= 2\Big(M[\tau](t,s)\,M[\tau](t{+}\delta,s) + M[-\tau'](t,s)\,M[-\tau'](t{+}\delta,s)\Big),
$$
$$
M[\text{leg}](a,b) = -\mathrm{Tr}\big[\Phi(a)\,\text{leg}(a,b)\,\Phi(b)\,\text{leg}(b,a)\big].
$$
The factor 2 is the two ways to pair the two source $\sigma$'s (both at $s$) with the two sink $\sigma$'s
(at $t,t{+}\delta$). Because $t,t{+}\delta\ne s$ for $dt\ge1$, EVERY leg is off-diagonal $\Rightarrow$ **no
equal-time contact, automatically vacuum-free** (as for PS $\langle\sigma^2 O_{2\sigma}\rangle_c=2M(t,s)M(t{+}\delta,s)$).
The overall $\tilde S$ minus is again innocuous ($M[-\tau']=M[+\tau']$: two legs per $M$, $(-1)^2=+1$), so
$M[-\tau'](a,b)=-\mathrm{Tr}[\Phi(a)\,\tau'(a,b)\,\Phi(b)\,\tau'(b,a)]$ with $\tau'=$ `tau_gw`.

Driver `fs_o2sigma_cross_claude.py`: per config build $M[\tau](a,b)$ and $M[\tau'](a,b)$ over the window,
form $S$-part $=\langle M[\tau]M[\tau]\rangle_s$, $\tilde S$-part $=\langle M[\tau']M[\tau']\rangle_s$, TOTAL
$=2(S{+}\tilde S)$; LINEAR, config jackknife. This is the GEVP off-diagonal $\langle\sigma_{FS}^2\,O_{\sigma\sigma}\rangle$.

### The cross has DISCONNECTED diagrams too (NM 2026-09-09) -- and for FS they do NOT vanish

The $2M\!\cdot\!M$ above is ONLY the two-meson (E-analog). The full $\langle O_{\sigma\sigma}(t,t{+}\delta)\,
\sigma^2(s)\rangle$ has the same contraction topologies as $\langle\sigma^2\sigma^2\rangle$ (4 $\sigma$
bilinears), but with the SINK time-split ($t,t{+}\delta$) and the source doubled at $s$. Enumerate them by a
brute-force Wick sum over the $4!=24$ pairings of the four vertices $p=(t,t{+}\delta,s,s)$: each permutation
$\pi$ factorizes into cycles, value $=(-1)^{\#\text{cycles}}\prod_{\text{cyc}}\mathrm{Tr}[\Phi(p_{c_0})\,
\text{leg}(p_{c_0},p_{c_1})\,\Phi(p_{c_1})\cdots\text{leg}(p_{c_k},p_{c_0})]$, with the equal-time contact
subtracted whenever a leg is coincident ($\text{leg}(a,a)$): $\tau(a,a)-\tfrac12 I$ for the $S$-part,
$-\tfrac12(\tau'+\tau)(a,a)$ for the $\tilde S$-part. FS $=$ $S$-part ($\tau$) $+$ $\tilde S$-part ($-\tau'$).

**Classification** (a permutation is CONNECTED iff some cycle contains both a sink vertex $\{t,t{+}\delta\}$
and a source vertex $\{s\}$ -- i.e. a fermion line bridges sink$\leftrightarrow$source):
- **CONNECTED (keep): A,B,C,D,E,G analogs.** E $=2M\!\cdot\!M$ (two mesons). A,B $=$ single connected loop
  through all 4. C,D,G carry a source/sink TADPOLE $D_S$ times a bridging loop/meson.
- **DISCONNECTED (drop): F,H,I,J analogs** -- sink piece $\times$ source piece, no bridging line.

In PS all the tadpole-carrying pieces (C,D,G and F,H,I,J) vanish because $\langle D_S^{PS}\rangle=0$
(contact-saturated), so the PS cross was cleanly $2M\!\cdot\!M$ (E only). **In FS the tadpole is live
($\langle D_S^{FS}\rangle=-0.96$), so C,D,G (connected, keep) and the disconnected F,H,I,J are all nonzero.**
The consistent treatment mirrors the diagonal: compute all, KEEP the connected (bridging) A,B,C,D,E,G, DROP the
disconnected F,H,I,J. Driver `fs_o2sigma_cross_diag_claude.py` (brute-force, connected vs disconnected).

## Plan once the contact is chosen
Write `fs_diag_corr_linear_claude.py`: per config compute `diags_pair(Phi, tau, s, t)` AND
`diags_pair(Phi, -taugw, s, t)` (the latter with the chosen FS contact on its diagonal), sum them per
diagram, translation-average, t-sum$\to$plateau subtract per jackknife, LINEAR per-diagram panels + TOTAL
$=\sum_i W_{10}(G_{10}[\tau]_i+G_{10}[-\tau']_i)$. Then the same B+E vs $2m_{PS}$ / A+C+D+G checks as PS.
