# Audit of `jj_corr_dilute_claude.cu` -- what it computes + consistency review

Critical read of the production diluted/master-field current-current driver against
`qed3int_v2-14.pdf` (source `qed3int_v2(7)/main.tex`). Companion design docs:
`jj_corr_dilute_impl_plan_claude.md`, `jj_local_stoch_estimator_design_claude.md`,
`jj_symmetries_masscases_claude.md`. All equation numbers below are **v2-14 PDF** numbers
(verified against the rendered PDF, not the older v2-11/12/13 numbering some code comments use).

File audited: `jj_corr_dilute_claude.cu` (942 lines).

---

## 0. One-paragraph summary

The program estimates, per gauge configuration and per stochastic hit, the **connected** and
(single-current trace pieces of the) **disconnected** vector and axial current-current correlators,
for three mass patterns (massless, flavor-breaking $m_F$, parity-breaking $m_P$), using one shared
$Z_2$ source $\eta$ and one shared forward solve $\varphi'=D_m^{-1}\eta$. It produces (i) the **exact
conserved-current** $K$ correlators (temporal `tp`, spatial `sp`), (ii) the **ultralocal** ($\sigma_a$)
current correlators `s1,s2,s3`, (iii) the **exact** and **local axial** $C_{A_{+-}}$, (iv) the
**disconnected single-current traces** $J(t)$, and (v) three **condensate** bilinear traces. Variance
reduction = spin $\times$ even/odd-time **dilution** (sum over disjoint patterns) + source-origin
**superposition** (master field), origins $t_0\in\{0,N_t/2\}$ folded downstream.

---

## 1. Definitions (operators, propagators, weights)

### 1.1 Overlap operators (code lines 292-294)
$$
D \equiv D_{\rm ov}\ (\text{massless}),\qquad
D_m \equiv D_{\rm ov}+m\ \ (\text{Eq. 3.59}),\qquad
\tilde D_{m_P}\equiv D_{\rm ov}+\frac{m_P}{1-m_P}\ \ (\text{Eq. 3.62}).
$$
`valence_mass = m = mass_re + i*mass_im` (line 234). Mass-pattern flags (lines 238-240):
`parity` $\equiv$ `|mass_im|>0 && |mass_re|~0` ($m_P\in i\mathbb R$); `flavor` $\equiv$
`|mass_re|>0 && |mass_im|~0` ($m_F\in\mathbb R$, now a dead flag -- see Issue I7).

### 1.2 Wick contractions / propagators (the "forward" and "backward" legs)
The two-point Wick rules of the PDF, by mass pattern:

| case | $\langle\xi\,\eta^\dagger\rangle$ (forward, $P$) | $\langle\eta\,\xi^\dagger\rangle$ (backward, $\bar P$) | eqs |
|---|---|---|---|
| massless | $D_{\rm ov}^{-1}$ | $D_{\rm ov}^{\dagger-1}$ | (3.50)/(3.51) |
| $m_F$ | $D_{m_F}^{-1}$ | $(D_{m_F}^{-1})^\dagger$ | (3.57)/(3.58) |
| $m_P$ | $D_{m_P}^{-1}$ | $(1+m_P)^{-1}(\tilde D_{m_P}^{-1})^\dagger$ | (3.60)/(3.61) |

In code, $P$ is realized as $\varphi'=D_m^{-1}\eta$ (line 542) and the source-leg solve
$D_m^{-\dagger}$ via the normal-equation `_ms` solver (estimator doc Sec. 7). The $m_P$ backward leg
uses $\tilde D_{m_P}^{-\dagger}$ (`phimm`, line 602). The kernel adjoint correction cancels the
$(1+m_P)^{-1}$ in the **exact** $V_{--}$/axial correlators (PDF text below Eq. 3.62), so the code does
**not** multiply $(1+m_P)^{-1}$ there; it **does** apply it to the **local** axial (line 899) and to
the condensate backward leg (lines 919-921), which carry no such cancellation.

### 1.3 Kernels and projection weights
Conserved-current kernel $K^{nn'}(t)=\partial D_{\rm ov}/\partial\theta_{nn'}(t)$ (`ConservedCurrent`,
line 297); temporal kernel $K^{t,t+1}(n)$ analogous. Weights (lines 358-369):
$$
w_{\rm tp}[n]=\frac{A_n}{\kappa_t(n)^2}\ (\text{Eq. 4.32}),\quad
w_{\rm sp}[\ell]=\frac{A_{nn'}}{\kappa^{(0)2}_{nn'}}\ (\text{Eq. 4.29}),\quad
w_{\rm site}[n]=A_n\ (\text{local; no }1/\kappa^2),
$$
with $A_n=$ `dual_areas`, $A_{nn'}=$ `link_volume`. Connected outputs are folded by $1/(4\pi)$
(`write_corr`, line 440); discs are RAW (`write_vec`, line 451).

---

## 2. The stochastic estimators (Wick $\to$ code)

Trace identity $\operatorname{tr}[M]=\mathbb E[\eta^\dagger M\eta]$. Source/sink asymmetry
(estimator doc Secs. 3-4): solve once on the **time-fixed source leg**, contract against the cheap
**shared sink** built on $\varphi'$.

**Exact vector $C_{V_{++},c}=\operatorname{tr}[K(t_0)\,P\,K(t)\,P]$** (Eq. 3.52). With
$P=D_m^{-1}$:
$$
\psi_{\rm tp}[n]=D_m^{-\dagger}\!\Big(\textstyle\sum_b K^\dagger(n,t_{0,b})\eta\Big),\qquad
\text{sink}=K(n,t)\,\varphi',\qquad
C_{\rm tp}(t)\mathrel{+}=w_{\rm tp}[n]\,\psi_{\rm tp}[n]^\dagger\,[K(n,t)\varphi'].
$$
Since $\psi^\dagger=\eta^\dagger K(t_0)D_m^{-1}$, the contraction $\psi^\dagger K(t)\varphi'
=\eta^\dagger K(t_0)D_m^{-1}K(t)D_m^{-1}\eta\xrightarrow{\mathbb E}\operatorname{tr}[K(t_0)PK(t)P]$.
**Verified correct** (lines 569-595). The summed-origin source makes the measured object
$C(t)=G(t)+G(t-N_t/2)$ (master field), folded in analysis.

**Exact disc $J(t)=\operatorname{tr}[D_m^{-1}K(t)]$** rides the same sink:
`Jtp[t] += w_tp[n]*eta.dag(kphi)` $=\eta^\dagger K(n,t)\varphi'$ (line 591). The disc correlator
$C_d(t_0,t)=J(t_0)J(t)$ (Eq. 3.53) is assembled downstream from the stored $J(t)$. **Correct.**

**Exact axial $C_{A_{+-}}=\operatorname{tr}[K(t_0)(1-D_{\rm ov}^\dagger)\bar P\,K^\dagger(t)
(1-D_{\rm ov})P]$** (Eq. 3.55 / 3.67). Code (lines 828-883): forward $\varphi'_A=D_m^{-1}\eta$ (line
833), $\chi=(1-D_{\rm ov})\varphi'_A$ (line 834, GW factor always massless); source
$\psi=X^{-\dagger}\,(1-D_{\rm ov})\,(\sum_b K^\dagger(t_{0,b})\eta)$ with $X=D_m$ (non-parity) or
$\tilde D_{m_P}$ (parity); sink $K^\dagger(t)\chi$. Gives $\bar P=D_m^{-\dagger}$ (massless/$m_F$) or
$\tilde D_{m_P}^{-\dagger}$ ($m_P$), $P=D_m^{-1}$. **Verified correct**, including the
massless-kernel-with-massive-propagators $m_F$ estimator (`jj_symmetries_masscases_claude.md` Sec. 3b)
and the $m_P$ form (Eq. 3.67). Only $C_{A_{+-}}$ is stored; $C_{A_{-+}}$ = reflection (Eq. 3.56),
reconstructed downstream.

**Local vector $\operatorname{tr}[\sigma_a(n,t_0)\,S\,\sigma_a(n,t)\,S]$, $S=D_m^{-1}$**
(lines 765-789). Localized source $\chi_{c,n}=D_m^{-\dagger}(\sigma_c\eta|_{(n,t_0)})$, local sink
$\sigma_c\varphi'$, diagonal-in-site dot product. **Verified correct** (estimator doc Sec. 4); both
legs forward $D_m^{-1}$. Local disc $J_c(t)=\operatorname{tr}[\sigma_c(n,t)D_m^{-1}]$ rides $\varphi'$.

**Local axial** (lines 799-820): same as local vector but the $t_0$ leg is the **backward**
propagator (`op_DmH` $\to\psi_A=D_m^{-1}\sigma_c\eta$, so $\psi_A^\dagger$ carries $D_m^{-\dagger}$);
i.e. $\operatorname{tr}[\sigma_a(n,t_0)\,D_m^{-\dagger}\,\sigma_a(n,t)\,D_m^{-1}]$ -- the local sibling
of Eq. 3.55 with the GW $(1-D_{\rm ov})$ "ripped off". For $m_P$ the $t_0$ leg is
$\tilde D_{m_P}^{-\dagger}$ and the $(1+m_P)^{-1}$ is applied at output (line 899). **Verified
consistent** with the design note.

**Condensates** (lines 904-928), area-weighted full traces:
$$
\sigma_{PS}=\langle\eta^\dagger\xi\rangle+\langle\xi^\dagger\eta\rangle\ (\text{Eq. 1.23}),\qquad
\sigma_{FS}=\langle\eta^\dagger\xi\rangle-\langle\xi^\dagger(1-D_{\rm ov}^\dagger)\eta\rangle\ (\text{Eq. 1.55}).
$$
with $\langle\eta^\dagger\xi\rangle=-\operatorname{tr}[A\,D_m^{-1}]$ (Grassmann sign, line 916), and
$\langle\xi^\dagger\eta\rangle$, $\langle\xi^\dagger(1-D_{\rm ov}^\dagger)\eta\rangle$ by mass case
(lines 917-925). **Verified correct**, including the subtle cyclicity+`[function of $D_{\rm ov}$
commute]` step that lets the massless/$m_F$ backward leg reuse $\varphi'$ via a conjugate.

---

## 3. Branch: configuration source -- `free_field` vs interacting (lines 459-468)

`free_field` $=$ `ens_dir.empty()`.
- **free**: one config $k=0$, $U=1$ (never read); `k_ckpoint=1,k_lo=0,k_hi=1`.
- **interacting**: loop `ckpoint_lat.k`, $k=k_{\min},k_{\min}+\texttt{stride},\dots<k_{\max}$;
  `U.read` per config; non-existent $k$ ends the loop (line 466).

No correctness issue. The output dir `esnid` encodes the mass and the ensemble basename (lines
389-396); dilution scheme is encoded in the subdir (`_s{0,1}td{td}`) so different schemes do not
collide.

---

## 4. Branch: dilution & superposition (source construction)

### 4.1 Dilution patterns (lines 533-537, `valence_claude.h:166-187`)
`n_spin_classes = spin_dilution ? NS : 1`; `time_dilution = td`; patterns $=$
`n_spin_classes * td`. `time_spin_dilution(rng,t_s,t_block=Nt/td,spin)` fills slices
$t\equiv t_s \pmod{td}$ for a single spin component; `time_dilution(...)` fills both spins. The
disjoint patterns tile the source, so the estimator is the **plain sum** over patterns (unbiased);
all hit-scope accumulators (lines 505-512) `+=` across patterns and are written **once after the
loop** (line 884+).

**Subtlety (benign, verified):** the source origins $t_0\in\{0,N_t/2\}$ are **both even**. For the
**local** current the source is localized exactly at $(n,t_0)$, so the odd time-dilution pattern
($t_s=1$) has zero local source and contributes nothing -- the full $t_0$-slice trace lives entirely
in the $t_s=0$ class, so summing patterns is still unbiased (not double-counted). For the **exact**
$K$ the temporal hop spans $\{t_0,t_0+1\}$ (opposite parities), so both classes contribute -- the
"standard 4-class partial reduction" caveat noted in the plan (Sec. 2). No bug.

### 4.2 Superposition / master field (lines 566-579, plan Sec. 3.5)
Source leg $=\sum_b K^\dagger(\cdot,t_{0,b})\eta$ (one solve per insertion); the measured correlator
is $C(t)=\sum_b G(t-t_{0,b})=G(t)+G(t-N_t/2)$, stored at **absolute** $t$ with `t0s={0,Nt/2}` for the
downstream fold. RNG seeded per `(config,hit)` from a string via `seed_seq` (lines 483-484),
reproducible from the stored `rng_seed`.

---

## 5. Branch: mass pattern (the big `if(parity)` fork)

The non-parity path (massless and $m_F$) is the validated production path. The `parity` ($m_P$) path
adds the $\tilde D_{m_P}$ "$--$" mirror and is **incomplete/broken under dilution** (see Issues).

### 5.1 Exact vector temporal `tp`
- **`(++)`** (lines 568-596, **all** mass cases): $C_{V_{++},c}$ (Eq. 3.52/3.63), superposed,
  pattern-summed into `Ctp`. Disc `Jtp` rides it. Written `h0/tp/Vpp` (line 887). Non-parity also
  writes `h0/tp/Vmm = conj(Vpp)` (Eq. 3.54 $C_{--}=C_{++}^*$).
- **`(--)`** (lines 599-653, **`if(parity)` only**): $C_{V_{--},c}=\operatorname{tr}[K^\dagger(t_0)
  \tilde D_{m_P}^{\dagger-1}K^\dagger(t)\tilde D_{m_P}^{\dagger-1}]$ (Eq. 3.65). **Logic verified
  correct**, but written per-origin `h0/t0_<b>/tp/Vmm` **inside** the dilution loop -- see Issues
  I1, I2.

### 5.2 Exact vector spatial `sp`
Same structure over links: `(++)` lines 655-685 (all cases), `(--)` lines 687-720 (parity only,
same per-origin-in-loop problem).

### 5.3 Exact axial `Apm`
Lines 828-883, all mass cases, superposed + pattern-summed, written `h0/axial/{tp,sp}/Apm` once after
the loop (lines 889-890). Parity uses $\tilde D_{m_P}$ sink leg. **Correct** (incl. parity).

### 5.4 Exact disc
`(++)` traces `Jtp/Jsp` ride the `(++)` sink (all cases), written `h0/disc/{tp,sp}/J` (lines
891-892). Parity `(--)` tilde traces `Jtil` (lines 726-757) -- written inside the loop, Issues
I1/I2.

### 5.5 Local vector `s1,s2,s3`
Lines 765-789, all cases; conn `h0/s{c}/Vpp`, non-parity `Vmm=conj(Vpp)`, disc `h0/disc/s{c}/J`
(lines 893-896). **Note:** for $m_P$ there is **no** local `(--)` channel computed -- Issue I3.

### 5.6 Local axial `s1,s2,s3`
Lines 799-820, all cases; `(1+m_P)^{-1}` applied for parity (line 899); written `h0/axial/s{c}/Apm`
(line 901). **Correct.**

### 5.7 Condensates
Lines 904-928 (all cases); `h0/condensate/{etadag_xi, xidag_eta, xidag_1mDdag_eta}`. **Correct.**

---

## 6. Issues found

### Severity: HIGH

**I1 -- Parity (`m_P`) crashes with more than one dilution pattern (duplicate HDF5 key).**
The parity `(--)` vector blocks (lines 647-651, 715-718) and the parity disc block (lines 753-755)
call `write_corr`/`write_vec` **inside** the dilution loop (`for t_s ... for sc`, lines 533-884).
HighFive `createDataSet` throws if the dataset already exists, so on the **second** pattern the key
`h0/t0_0/tp/Vmm` (etc.) already exists and the program aborts. The default invocation (no spin
dilution, `--time-dilution 2`) already has 2 patterns, and the spin-only recommendation
(`--spin-dilution --time-dilution 1`) has 2 spin patterns -- so **every** realistic $m_P$ run
crashes. Because the abort happens mid-hit, the per-hit `.tmp` file never gets the `complete`
sentinel (good -- no corrupt "complete" file), but no $m_P$ data is produced. The condensate and the
exact axial for $m_P$ are also lost since the crash precedes the post-loop writes.
*Fix direction:* convert the `(--)` and parity-disc blocks to the same superposed + hit-scope
accumulator + write-once-after-loop pattern as the non-parity path (the plan already flags this as
"PARITY-TODO", impl plan Sec. 9).

**I2 -- Parity output is structurally inconsistent with the non-parity output even at 1 pattern.**
Even if I1 is sidestepped (single pattern), $m_P$ writes `Vpp` superposed at `h0/tp/Vpp` (one folded
$C(t)=G(t)+G(t-N_t/2)$) but `Vmm` per-origin at `h0/t0_0/tp/Vmm`, `h0/t0_1/tp/Vmm` (un-superposed).
The two channels are then in different representations and cannot be combined; for $m_P$ the physics
of interest is exactly $C_{V_{++}}\neq C_{V_{--}}^*$ (the parity signal,
`jj_symmetries_masscases_claude.md` Sec. 3c), which requires `Vpp` and `Vmm` in the **same** frame.

### Severity: MEDIUM

**I3 -- Local vector `(--)` missing for `m_P`.** The local vector (lines 765-789) only builds the
`Vpp` channel ($D_m=D_{m_P}$ both legs); there is no $\tilde D_{m_P}$ local `(--)`. The write guard
`if(!parity)` (line 895) silently drops local `Vmm` for $m_P$. Whether this matters depends on
whether the local channel is used for $m_P$ at all (it is approximate/UV by design), but it should be
a documented omission, not a silent one.

**I4 -- `Ctp`/`Csp` name shadowing in the parity blocks.** Lines 631 and 703 declare
`std::vector<std::vector<Complex>> Ctp(...)` / `Csp(...)` that **shadow** the hit-scope accumulators
`Ctp`/`Csp` (lines 505). It is currently harmless (the outer `(++)` accumulation at lines 592/681
runs before the inner scope, and the post-loop write uses the outer one), but it is a footgun: a
future edit touching the inner blocks could read/write the wrong `Ctp`. Rename the local parity
buffers (e.g. `Ctp_mm`).

### Severity: LOW / cosmetic

**I5 -- Equation-number drift in code comments.** Lines 290-291 cite "Eq. 3.60" for $D_m=D_{\rm ov}+m$
(actual: **3.59**) and "Eq. 3.63" for $\tilde D_{m_P}=D_{\rm ov}+m_P/(1-m_P)$ (actual: **3.62**).
These look like leftovers from the v2-12/13 numbering. The estimator-side citations (3.50/3.51,
3.55, 3.57/3.58, 3.60/3.61, 3.65/3.67, 3.68, 4.29, 4.32, 1.23, 1.55) are all **correct** for v2-14.

**I6 -- Duplicate LaTeX label in the source paper (not the code).**
`qed3int_v2(7)/main.tex` lines 2463 and 2498 both `\label{eq:est_spatially_projected_JJ}`, so any
`\eqref` to it is ambiguous (resolves to the second). Worth fixing in the manuscript since the audit
relies on these numbers.

**I7 -- Dead `flavor` flag.** `flavor` is computed (line 240) then `(void)flavor` (line 241); the
$m_F$ axial now uses the $D_m$ path (same as massless), with the old massless-$D_{\rm ov}$ branches
left commented (lines 831-832, 842-848, 864-868). This is intentional (commented A/B per house
style) but means the `flavor` branch reads as live when it is not -- a one-line comment at the flag
would help.

**I8 -- `TOL_OUTER=1e-5` (line 90).** Loosened from the older `1e-8`; the comment justifies it for
the small-norm sink RHS. Not an error, but it is a global knob shared by **all** solves (condensate,
axial, local) -- worth a sanity check that the axial GW-dressed solve tolerates it (the axial RHS
$(1-D_{\rm ov})\rho$ is not small-norm in the same way as the bare sink). Flagging for awareness, not
as a known bug.

---

## 7. What is verified correct (so it is on record)

- Non-parity (massless, $m_F$) vector `tp`/`sp` `(++)`/`(--)=conj`, exact axial `Apm`, exact disc,
  local `s1,s2,s3` vector + axial, and all three condensates: estimator algebra, propagator legs,
  Grassmann signs, weights, and $1/(4\pi)$ folding all match the PDF (Sec. 1-2 above) and the
  deterministic references cited in `MEMORY.md`.
- Dilution summation is unbiased (Sec. 4.1), including the benign odd-class-vanishes behaviour for
  the local current.
- Superposition / master-field source and the `C(t)=G(t)+G(t-N_t/2)` storage (Sec. 4.2).
- Parity `(++)` vector, parity exact axial, parity condensate **algebra** is correct (Eq. 3.63/3.67);
  only the *output plumbing* for the parity `(--)` and parity-disc blocks is broken (I1/I2).

**Next:** the actionable item is I1/I2 -- refactor the three parity blocks (lines 599-653, 687-720,
726-757) to superposed hit-scope accumulators written once after the dilution loop, mirroring the
non-parity path, so $m_P$ runs at all and its `Vpp`/`Vmm` live in the same frame.

---

## 8. Resolution (I1/I2) (2026-06-13) -- plan `mP_parity_refactor_impl_plan_claude.md`

### 8.1 Files + line ranges touched (`jj_corr_dilute_claude.cu`)
| block | new lines | what changed |
|---|---|---|
| parity hit-scope accumulators | `:517-519` | NEW: `Ctp_mm`,`Csp_mm` (exact-K V--), `JtpT_til`,`JspT_til` (disc), `Cs_mm[3]` (local V--) |
| `(--)` vector tp | `:606-640` | per-origin write-in-loop -> SUPERPOSED summed source + accumulate `Ctp_mm`, NO write; condensate (`:611-616`) kept; inert ylm dropped |
| `(--)` vector sp | `:675-699` | same: summed source + accumulate `Csp_mm`, NO write |
| `(--)` disc tilde | `:706-730` | local `JtpT/JspT` (written-in-loop) -> hit-scope `JtpT_til/JspT_til` `+=`, NO write |
| local vector V-- (NEW, I3) | `:761-775` (in the local block `:733-780`) | parity branch: `psi_mm=tilde D^{-1}(sigma_c eta)`, sink `sigma_c phimm` -> `Cs_mm[c-1]` |
| post-loop writes | `:880-890` (exact V--/Jtil), `:893-897` (local V--) | the ONLY parity writes; after the dilution loop |
| cosmetic I5 / I7 | `:289-290` / `:241-242` | eq# 3.59/3.62; dead-`flavor` note |
Non-parity (massless/$m_F$) path: BYTE-UNCHANGED. Old in-loop write lines are kept commented beneath the
new post-loop writes (`:885-888`) for A/B; the NEW (audit I3) local-V-- write is marked at `:896`.

### 8.2 New key layout for the `(--)` and parity-disc outputs -- CLOSES I2
All written ONCE after the dilution loop at the **`h0/` level (NO `t0_<b>` nesting)**, the SAME superposed
frame ($C(t)=G(t)+G(t-N_t/2)$, absolute $t$) as `Vpp`:
- `h0/tp/Vmm`, `h0/sp/Vmm`        (exact-K V--, Eq. 3.65; tilde D, NO $(1+m_P)$ factor -- kernel-cancelled)
- `h0/disc/tp/Jtil`, `h0/disc/sp/Jtil`   (tilde disc traces)
- `h0/s{1,2,3}/Vmm`               (local V--, $(1+m_P)^{-2}$, Eq. 3.65 local)
- (already same-frame: `h0/tp/Vpp`,`h0/sp/Vpp`,`h0/s{c}/Vpp`,`h0/axial/...`, condensate)
The OLD per-origin `h0/t0_<b>/{tp,sp}/Vmm` keys are GONE. `Vpp`/`Vmm` now combinable for the parity signal
$C_{V++}\neq C_{V--}^*$.

### 8.3 Accumulator confirmation -- CLOSES I1
- exact-K `(--)`: `Ctp_mm[t] += w_tp[n]*psi.dag(kphi)` (`:638`), `Csp_mm[t] += ...` (`:697`) -- hit-scope,
  `+=` across ALL dilution patterns, written once (`:881-882`).
- parity disc: `JtpT_til[t] += ...` (`:717`), `JspT_til[t] += ...` (`:727`) -- hit-scope, `+=` across
  patterns, written once (`:883-884`).
- local `(--)`: `Cs_mm[c-1][t] += ...` (`:774`) -- hit-scope, written once (`:896`).
No `write_corr`/`write_vec` remains inside the `for t_s ... for sc` loop (`:540-875`) -> no duplicate-key
abort with >=2 patterns. (I4: the shadowing local `Ctp`/`Csp` buffers are removed; renamed to `Ctp_mm`/`Csp_mm`.)

### 8.4 Deviations from "Fix direction" + I3 status
- **No deviation** on I1/I2: did exactly the superposed + hit-scope + write-once refactor.
- **Bonus fix:** the old per-origin `(--)` also had a latent `psi_tp[ITP(n,b)]=psi_tp[n]` aliasing bug
  (`ITP` ignores `b`, so it kept only the last origin); the superposed rewrite removes it.
- **I3 -- ADDRESSED (computed, not left).** Local V-- for $m_P$ is now built (tilde both legs, $(1+m_P)^{-2}$).
  Note the differing $(1+m_P)$ powers (bare currents, no kernel cancellation): local V-- $(1+m_P)^{-2}$,
  local axial $(1+m_P)^{-1}$, condensate backward legs $(1+m_P)^{-1}$; exact-K V--/axial carry NONE
  (kernel-cancelled). Local axial $m_P$ was ALSO fixed to tilde (dilute + `jj_local_axial_deter`), and the
  exact-K determ got tilde-D (`jj_exact_diag_deter_free` V--, `jj_exact_axial_deter_free`).

### 8.5 Validation note (correction to the reviewer's)
`parity` is **mass-driven** (`mass_im != 0 && mass_re ~ 0`, `:238`), INDEPENDENT of the gauge field -- there
is NO `free_field` guard on the `(--)` blocks. So a **FREE run with `--mass-im`** DOES exercise the fix
(parity=true, tilde solves run; and tilde D $\neq D_m$ even at $U=1$, so the signal is nonzero). Thus the
crash/frame fix is provable WITHOUT an interacting run:
```
nvcc <flags> jj_corr_dilute_claude.cu -o jj_corr_dilute_mPtest.o
./jj_corr_dilute_mPtest.o --nhits 1 --n-t0 2 --spin-dilution --time-dilution 2 --mass-im 0.1 2>&1 | tee mPtest_claude.log
```
This is 2 patterns (exactly what used to crash). Proof of fix = it completes and the per-hit `.h5` contains
`h0/{tp,sp}/Vmm`, `h0/disc/{tp,sp}/Jtil`, `h0/s{c}/Vmm` at the `h0/` level (no `t0_<b>` Vmm). The production
`run_jj_dilute_valid_mP_claude.sh` is itself a free `--mass-im 0.1` run, so it doubles as the proof; teeing
its log lets you read back (a) no duplicate-key abort and (b) the key layout.
