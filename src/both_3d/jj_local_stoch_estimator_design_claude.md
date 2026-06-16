# Design of the stochastic LOCAL-current estimators (diagonal scheme)

Companion to `jj_stoch_local_impl_plan_claude.md`. This explains *why* the estimators are built the way
they are, from the stochastic trace identity down to the per-site solve structure. Scheme: **diagonal**
(user-confirmed 2026-06-12); the deterministic reference is `jj_local_deter_claude.cu`
(`corr_deter_local_L<L>`).

---

## 1. What we are estimating

The deterministic local connected correlator is an explicit trace over the whole lattice:

$$
G_a(t)=\sum_n w_n\,\mathrm{tr}\big[\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S\big],
\qquad S\equiv D_m^{-1},\quad w_n=\texttt{dual\_areas}[n].
$$

$\sigma_a(n,t)$ denotes the Pauli matrix $\sigma_a$ inserted at site $n$, timeslice $t$ -- i.e. $\sigma_a$
times the projector onto that single site. The defining feature of this observable (the "diagonal sum",
Eq. 4.29) is that **the same site $n$ appears at both $t_0$ and $t$**. There is no second, independent
site sum; that would be the $q=0$ double sum, which is a different (and here unwanted) object.

Channels: $a=1,2$ are spatial ($G_s=s_1+s_2$, Eq. 4.28), $a=3$ is temporal ($G_t$, Eq. 4.31).

---

## 2. The stochastic trace identity

A single random source $\eta$ with $\mathbb{E}[\eta_i\,\eta_j^*]=\delta_{ij}$ gives the unbiased estimator

$$
\mathrm{tr}[M]=\mathbb{E}\big[\eta^\dagger M\,\eta\big].
$$

Applied per site $n$,

$$
\mathrm{tr}\big[\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S\big]
=\mathbb{E}\big[\eta^\dagger\,\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S\,\eta\big].
$$

Everything below is about evaluating the right-hand side with as few linear solves as possible, while
keeping it pinned to the diagonal in $n$.

---

## 3. Which side carries the solve -- the source/sink asymmetry

Read the operator string $\eta^\dagger\,\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S\,\eta$ from the right:

- **Right end:** $S\eta=\texttt{phi'}$ -- a single solve, **shared by everything** (all $n$, all channels,
  connected and disconnected).
- Continuing from the right would hit $\sigma_a(n,t)\,\texttt{phi'}$ and then another $S[\cdots]$;
  applying $S$ to a source localized at the *sink* $(n,t)$ would cost **one solve per sink** $(n,t)$. Bad.
- **So fold from the left instead.** Define

$$
\texttt{chi}_{n,a}\equiv S^\dagger\,\sigma_a(n,t_0)\,\eta
=D_m^{-\dagger}\big(\sigma_a\,\eta|_{(n,t_0)}\big),
\qquad
\texttt{chi}_{n,a}^\dagger=\eta^\dagger\,\sigma_a(n,t_0)\,S .
$$

The source $\sigma_a\,\eta|_{(n,t_0)}$ lives at a **single point** $(n,t_0)$ and is **independent of the
sink time** $t$ -- one solve per $(n,t_0)$, reusable for all $t$.

The estimator then collapses to a **local** contraction:

$$
\eta^\dagger\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S\,\eta
=\texttt{chi}_{n,a}^\dagger\big[\sigma_a(n,t)\,\texttt{phi'}\big]
=\texttt{chi}_{n,a}(t,n)^\dagger\,\sigma_a\,\texttt{phi'}(t,n),
$$

a 2-spinor dot product at $(n,t)$. The sink costs **no solve** because $\sigma$ is ultralocal. This
asymmetry -- solve on the time-fixed source, free local sink -- is the core efficiency idea (same as the
conserved-current $K$ code).

---

## 4. Why localizing the source gives the diagonal (the crux)

The two $\eta$'s play different roles:

- the **left** $\eta$ enters only through $\sigma_a\,\eta|_{(n,t_0)}$ -- it uses $\eta$ **only at $(n,t_0)$**;
- the **right** $\eta$ is inside $\texttt{phi'}=S\eta$ -- it uses $\eta$ **everywhere**.

Taking the expectation, the contraction $\mathbb{E}\big[\eta|_{(n,t_0)}\,\eta^\dagger\big]$ is the
**projector onto $(n,t_0)$**. It forces the right $\eta$ to sit at the *same* point $(n,t_0)$, pinning both
propagator legs to connect $(n,t_0)\leftrightarrow(n,t)$:

$$
\mathbb{E}\big[\texttt{chi}_{n,a}^\dagger\,\sigma_a(n,t)\,\texttt{phi'}\big]
=\mathrm{tr}\big[\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S\big].
$$

Any contribution where the right $\eta$ sat at a different site $n'\neq n$ needs
$\mathbb{E}[\eta_{(n,t_0)}\,\eta_{(n',t_0)}^*]=0$: it vanishes in expectation and only contributes
**variance**. So *localizing the source automatically projects onto the diagonal* -- no explicit
restriction of the contraction is required. This is also why the per-site solves are unavoidable: each
$n$ needs its own localized source.

**Contrast (the meson wall).** If the source were $\sigma_a\,\eta|_{\text{wall }t_0}$ (all sites on the
$t_0$ slice), then $\mathbb{E}[\eta_{\text{wall}}\,\eta^\dagger]$ projects onto the whole slice, the right
$\eta$ can land on any $n'$, and one gets the **double sum** $\sum_{n,n'}$ ($q=0$ correlator):

$$
\text{local/point source}\ \Longleftrightarrow\ \text{diagonal }\sum_n,
\qquad
\text{wall source}\ \Longleftrightarrow\ \text{double sum }\sum_{n,n'} .
$$

This is exactly the A-vs-B fork, now read off from the contraction algebra. We want A (diagonal).

---

## 5. The ylm tower is intentionally extended

The ylm insertion is itself a site sum, $\Sigma_{lm}=\sum_n A_n\,Y_{lm}(\hat n)\,\sigma_3(n)$ with
$A_n=\texttt{dual\_areas}[n]$, and the observable

$$
g_l(t)=\frac{1}{2l+1}\sum_{m}\mathrm{tr}\big[\Sigma_{lm}(t_0)\,S\,\Sigma_{lm}(t)\,S\big]
$$

genuinely wants the full $\sum_{n,n'}A_nA_{n'}Y_{lm}(\hat n)Y_{lm}(\hat n')$ structure. So here a
**wall-type source is correct**:

$$
\Sigma_{lm}(t_0)\,\eta=\sum_n A_n Y_{lm}(\hat n)\,\sigma_3\,\eta|_{(n,t_0)}
=\big(\eta\text{ on the }t_0\text{ wall}\big)\xrightarrow{\ \sigma_3\ }\xrightarrow{\ \times A_nY_{lm}\ }.
$$

Conveniently `mult_Ylm_real(l,m,base)` already multiplies each site by $\texttt{dual\_areas}\cdot Y_{lm}
= A_n Y_{lm}$, so the source is literally `eta_wall.mult_sigma3().mult_Ylm_real(l,m,base)`. One solve per
$(l,m,t_0)$; $\sum_{l=0}^{2}(2l+1)=9$ solves per $t_0$. The sink mirrors it:
$\texttt{sig3phi\_ylm}=(\sigma_3\,\texttt{phi'})$ then `.mult_Ylm_real(l,m,base)`, contracted per sink slice.

---

## 6. Disconnected rides on `phi'` alone

$$
J_a(t)=\sum_n w_n\,\mathrm{tr}[\sigma_a(n,t)\,S]
=\mathbb{E}\Big[\sum_n w_n\,\eta(t,n)^\dagger\,\sigma_a\,\texttt{phi'}(t,n)\Big].
$$

Only the shared $\texttt{phi'}$ is needed -- **no source solves**. That is why the disconnected piece is
the cheap, many-config path and lives in its own file. The ylm disc is the same with the $A_nY_{lm}$
weight on $\sigma_3$.

---

## 7. Realizing $S$ and $S^\dagger$, and the sign

Overlap inverses via the normal equations (multishift `_ms` for the squared operator
$D_m^\dagger D_m$):

$$
S\,\eta=D_m^{-1}\eta:\quad (D_m^\dagger D_m)\,z=D_m^\dagger\eta\ \Rightarrow\ z=D_m^{-1}\eta,
$$
$$
S^\dagger x=D_m^{-\dagger}x:\quad (D_m^\dagger D_m)\,z=D_m\,x\ \Rightarrow\ z=D_m^{-\dagger}x.
$$

In code: `op_DmH.from_cpu(tmp,eta)` then `op_Dmsq.solve(phi',tmp)` for $S\eta$;
`op_Dm.from_cpu(tmp,x)` then `op_Dmsq.solve(chi,tmp)` for $S^\dagger x$.

Because the operator ordering is kept identical to the deterministic trace and the same $P=D_m^{-1}$ is
used, the **sign and normalization match exactly**: the stochastic `Vpp` should land on the deterministic
`corr_deter_local_L1` within the jackknife error, with no fudge factor. (As in the determ code,
`write_corr` folds $1/(4\pi)$ and `Vmm=conj(Vpp)` for the non-parity case.)

---

## 8. Variance and resampling

One source $\eta$ per hit; the off-diagonal $\eta$-noise (the $n'\neq n$ terms in Sec. 4) averages out
over hits. Errors are jackknife over **hits** in the free field (one gauge config), and over **configs**
in the interacting theory -- the same convention as the $K$ stochastic notebooks.

---

## 9. Cost summary (per hit)

$$
\underbrace{1}_{\texttt{phi'}}
\;+\;\underbrace{3\,n_{\text{sites}}\,n_{t_0}}_{\text{conn }s_1,s_2,s_3}
\;+\;\underbrace{9\,n_{t_0}}_{\text{conn ylm}}
\quad\text{solves (connected file)};
\qquad
\underbrace{1}_{\texttt{phi'}}\ \text{solve (disconnected file).}
$$

At $L=1$ ($n_{\text{sites}}=12$, $n_{t_0}=2$): $\approx 1+72+18=91$ connected solves per hit. The factor-3
over channels and the per-site loop can later be batched with the mrhs block solver (as in
`jj_corr_block_t`), but v1 keeps the simple per-source loop for clarity.
