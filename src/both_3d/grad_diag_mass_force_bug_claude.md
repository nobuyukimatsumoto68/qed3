# Diagonal-mass HMC force bug: the $(1+M^*)$ fold is wrong for site-varying $m_L$ ($L>1$)

**Status (2026-06-19): FIXED + FD-VALIDATED.** L=2 massive force-vs-FD dropped ~1e-4 -> ~5e-9..3e-8
(real m=0.1 AND imag m=0.1i), matching massless; the grad VALUE corrected (real-m sp2 3.764e-2 ->
3.7536e-2 = fd) and L2_default == L2_gradl4 byte-identical. Operator checks unchanged. REMAINING:
(a) re-run HMC test at normal tmax -> dH ratio should now be ~4 (the tmax_sc=0.2 workaround can be reverted);
(b) remove the now-unused d_MA/d_MB/d_dotAM/d_dotBM buffers; (c) update memory. Detail of the fix below.

---
(historical) FIX IMPLEMENTED (Sec 5'): `precalc_grad`
(default + `_ms`) now does the extra multishift solve `W_m = R_m X^dag(M_mass eta)`, overwrites
`d_Ys[m] -> hat Y_m = Y_m + mass_coeff*W_m`, and sets `d_eta_bra = (1+m_L)eta`; `grad`/`l1`/`l2`/`l4`
dropped the `conj(mass_coeff)*M_mass` fold and use `hat Y` (bra) + `d_eta_bra` (m=0). Massless = identical.
New members `d_eta_bra/d_Mmeta/d_Ws`. VALIDATION: rebuild+run `tmp_claude.sh`; L=2 force-vs-FD massive
rows must drop `~1e-4 -> ~1e-8` (default AND grad_l4); then HMC dH-scaling ratio `~4` at normal tmax.
(Old fold buffers `d_MA/d_MB/d_dotAM/d_dotBM` now unused -- remove after FD confirms.)
Scope: the HMC **force** (`grad_deviceAsyncLaunch` + `grad_l1/l2/l4` + `precalc_grad`) in
`includes/overlap_wmass_claude.h`. The **operator** (`mult`/`adj`/`DHD` + `_ms`, and
`BlockedMat`) is verified correct (adjointness ~1e-10, BlockedMat-vs-MatPoly = 0.0, obsolete-vs-prod
1e-16) — this is force-only. Production HMC uses `#define GRAD_L4`.

Context: measure-weighted diagonal mass $m_L=\text{diag}(m\,A_y/\bar a_s)=\text{mass\_coeff}\cdot
M_\text{mass}$, $M_\text{mass}=\text{diag}(A_y/\bar A)$ (the volume factor $V$, real diagonal),
$\text{mass\_coeff}=m\,\bar A/\bar a_s$. See `mass_measure_factor_impl_plan_claude.md`.

## 1. Symptom

- HMC `test_hmc_diag_mass_claude.cu` **Test 2 (dH ~ $\tau^2$ scaling)** FAILS for $m\neq0$ at $L=2$:
  ratio $\approx 1.3$ ($m=0.1$) / $2.4$ ($m=0.1i$) vs the expected $\approx 4$; **massless passes**
  (ratio $3.90$). Reversibility ($\sim$1e-10) and acceptance ($<\!e^{-dH}\!>\approx1$) pass.
- A constant-ish $dH$ that does NOT shrink $\sim\!4\times$ when $n_\text{steps}$ doubles, and that
  PERSISTS as $\tau\to0$, is the signature of a **non-conservative force**: force $\neq dS/dU$. The
  integrator conserves the shadow Hamiltonian of whatever force it uses; if that force is not
  $dS/dU$, the true $H$ drifts by a $\tau$-independent amount $\Rightarrow$ a $dH$ floor.

## 2. Evidence it is the force (not the operator, not a grad variant)

From the operator test `test_diag_mass_l1_claude.cu` force-vs-finite-difference
(`grad` vs central difference of $S=\phi^\dagger(D_m^\dagger D_m)^{-1}\phi$):

| L | massless $|grad-fd|$ | massive ($m=0.1$) $|grad-fd|$ |
|---|---|---|
| 1 | ~6e-8 (clean, link-uniform) | ~6e-8 (clean) |
| 2 | ~1e-9..1e-5 | **~1e-4..1.7e-4, link-DEPENDENT (systematic)** |

- The $L=1$ massive gap is FD-precision-limited (~6e-8) $\Rightarrow$ force is correct at $L=1$.
- The $L=2$ massive gap is `~1e-4`, **link-dependent and systematic** (not the random size of FD
  noise) $\Rightarrow$ a real force error at $L>1$. (Earlier mis-attributed to "FD noise" — wrong.)
- `grad_l4` $\equiv$ reference `grad` to machine precision (L=1 obsolete-vs-prod `2.2e-16`; L=2
  `grad` values byte-identical between the `L2_default` and `L2_gradl4` builds). So this is NOT a
  block-variant bug — ALL grad variants share the same (incorrect) mass fold.

## 3. Root cause

The correct pseudo-fermion force, with GW $D_\text{ov}^\dagger D_\text{ov}=D_\text{ov}+D_\text{ov}^\dagger$,
$dD_m/dU = dD_\text{ov}/dU \equiv K$ (the conserved current, **mass-independent**), and $dM/dU=0$:
$$
\frac{dD_m^\dagger D_m}{dU} = (1+M^*)K + K^\dagger(1+M),
\qquad
\frac{dS}{dU} = -\eta^\dagger\big[(1+M^*)K + K^\dagger(1+M)\big]\eta
= -2\,\mathrm{Re}\,\eta^\dagger(1+M^*)K\,\eta ,
$$
with $\eta=(D_m^\dagger D_m)^{-1}\phi$ (the existing massive solve). So $(1+M^*)$ sits **immediately
after the outer bra $\eta^\dagger$** — nowhere else. Equivalently $-2\,\mathrm{Re}\langle(1+M)\eta\,|\,K\eta\rangle$.

**What the code does instead (the bug).** The massless grad computes $\langle\eta|K\eta\rangle$ in a
factored Zolotarev-pole form, with $X=M_\text{DW}/\lambda_\text{max}$, $R_m=(XX^\dagger-k^2/cp[m])^{-1}$,
$Z_m=R_m\eta$ (ket), $Y_m=R_m X^\dagger\eta$ (bra), and a per-link COO $W^{wz}$ ("coo"). The pole
inner products are e.g. $\langle\text{coo}\,Y_m\,|\,X Z_m\rangle$. The current mass fold multiplies
$$
\langle\text{coo}\,Y_m\,|\,M_\text{mass}\,(X Z_m)\rangle ,
$$
i.e. it inserts $M_\text{mass}$ on $X Z_m$, **deep inside** the contraction (after $X$, after $R_m$).
Because $V=M_\text{mass}$ does NOT commute with $X$ or $R_m$, this is the WRONG slot. It only
coincides with the external-bra $(1+M^*)$ at $L=1$ where $M_\text{mass}=\mathbb 1$ (commutes with
everything) — which is exactly why $L=1$ is exact and $L>1$ breaks. It also reduces to the validated
SCALAR $(1+m^*)$ fold at $L=1$.

**The $m=0$ (Term-A/explicit-bra) term is already correct:** there the bra $\eta$ is explicit, so
`<eta | M_mass coo(...)>` $=\eta^\dagger M_\text{mass}(\dots)=(M_\text{mass}\eta)^\dagger(\dots)$
(single Hermitian insertion). Only the **pole terms** ($m>0$) are wrong.

## 3b. Grounding: Eq.(D.7) of arXiv:1606.04109 (Karthik-Narayanan) + the generalization

The code's massless force IS Eq.(D.7) of arXiv:1606.04109 App.D (HMC force, $N_f=1$ massless
overlap). With $\Psi\equiv(C_o^\dagger C_o)^{-1}\phi$ (= our `eta`), $X'=\partial X/\partial\theta$
(= our `coo`/$W^{wz}$), $R_k=(X^\dagger X+p_k)^{-1}$ (Zolotarev resolvent, multishift), and
$$Y_k=R_kX^\dagger\Psi,\quad Z_k=r_kR_k\Psi,\quad \tilde Y_k=XY_k,\quad \tilde Z_k=XZ_k,$$
$$F_\mu=-\tfrac12\mathrm{Re}\Big\{\Psi^\dagger X'\big(\textstyle\sum_k Z_k\big)-\sum_k\big(\tilde Y_k^\dagger X'Z_k+\tilde Z_k^\dagger X'Y_k\big)\Big\}.\quad\text{(D.7)}$$
Symbol map to `overlap_wmass_claude.h`: $\Psi\!\to$`d_eta`; $Y_k\!\to$`d_Ys[m]`; $Z_k\!\to$`d_Zs[m]`
($r_k\!\to$`A[m]`); $\tilde Y_k=XY_k\!\to$`d_XYs[m]` (`X Y_m`); $\tilde Z_k=XZ_k\!\to$`d_XZs[m]`
(`X Z_m`); $X'\!\to$`coo`. The three D.7 pieces $=$ {the $m{=}0$ term, term b $\langle XY_m|\mathrm{coo}\,Z_m\rangle$,
term a $\langle\mathrm{coo}\,Y_m|XZ_m\rangle$}. The **tilde vectors $\tilde Y_k,\tilde Z_k$ are the
BRA-side vectors** (they carry $\Psi^\dagger$); $Y_k,Z_k$ are ket-side.

**Generalization to our setting (NM):** the paper is massless; our $C_o=D_m=D_\text{ov}+M$,
$M=\text{mass\_coeff}\cdot V$, $V=M_\text{mass}$ (diagonal, "special setting"). The mass enters
$F=\mathrm{Re}[(C_o\Psi)^\dagger\dots]$ via the OUTER bra $C_o\Psi=(1+M)\,$(massless bra structure),
i.e. by pre-multiplying $(1+M)$ onto the **bra-side tilde copies** (and the explicit $\Psi$):
$$
F_\mu^{(m)}=-\tfrac12\mathrm{Re}\Big\{\big((1{+}M)\Psi\big)^\dagger X'\big(\textstyle\sum_k Z_k\big)
-\sum_k\big(((1{+}M)\tilde Y_k)^\dagger X'Z_k+((1{+}M)\tilde Z_k)^\dagger X'Y_k\big)\Big\}.
$$
So: make NEW copies $(1{+}M)\tilde Y_k$, $(1{+}M)\tilde Z_k$, $(1{+}M)\Psi$ — a cheap **diagonal
multiply on the already-solved, X-applied vectors** (`X Y_m`, `X Z_m`, `eta`) — **no new solve**.
For SCALAR $m$ ($M=m\mathbb 1$) this is exactly the old validated $(1+m^*)$ fold (the scalar commutes
through everything). For diagonal $V$ it does NOT (the flagged non-commutativity), and the placement
**on the tilde (post-$X$) bra vectors** is what the current code gets wrong (see Sec 5).

## 4. The constraint (NM): NO additional solve

Naively moving $M_\text{mass}$ to the external bra of pole term a,
$\eta^\dagger\,X R_m\,\text{coo}^\dagger\,X Z_m \to (M_\text{mass}\eta)^\dagger X R_m\,\text{coo}^\dagger X Z_m
= \langle \text{coo}\,R_m X^\dagger(M_\text{mass}\eta)\,|\,X Z_m\rangle$,
*appears* to need $R_m$ on a NEW right-hand side $X^\dagger(M_\text{mass}\eta)$ — a third multishift
solve in `precalc_grad`. **NM asserts no additional solve is needed**: the exact force is
reproducible from the EXISTING $\eta$ (and existing $Z_m=R_m\eta$, $Y_m=R_m X^\dagger\eta$); it is a
**multiplication-ordering** issue, non-commutative in $V$. So $M_\text{mass}$ must go in a slot that
uses only already-computed vectors + cheap (diagonal/matvec) ops.

## 5'. CONFIRMED FIX (NM 2026-06-19): one extra multishift solve on the bra

FD data settled it: massless L=2 force-vs-FD is CLEAN (~5e-9..3e-8); BOTH real `m=0.1` (~5e-6..1.7e-4)
and imag `m=0.1i` (~3e-5..1.6e-4) are SYSTEMATICALLY wrong. The `(1+M)`-on-tilde (Sec 5) is provably
a NO-OP for real $m$ ($M_\text{mass}$ Hermitian => moving it across a single dot changes nothing), so
it cannot be the fix. NM confirmed: the additive operator needs the bra `(1+M)eta` carried THROUGH the
resolvent -> **one extra multishift solve** (massless path keeps the massless multishift).

Exact force ($D_m=D_\text{ov}+M$, $M=m_L=\text{mass\_coeff}\cdot M_\text{mass}$, $K$ mass-indep, GW):
$$F=\eta^\dagger(K^\dagger D_m+D_m^\dagger K)\eta=2\mathrm{Re}\,[(D_m\eta)^\dagger K\eta]
=2\mathrm{Re}\,[((D_\text{ov}+M)\eta)^\dagger K\eta].$$
The massless D.7 (Sec 3b) is the $(D_\text{ov}\eta)$ part with the massive $\eta$. The extra
$(M\eta)$ part has bra$\to M\eta$, which enters the bra-built vectors $Y_k=R_kX^\dagger\eta$ and the
explicit $\Psi^\dagger$ (the ket-built $Z_k,\tilde Z_k$ are UNCHANGED). Net: replace the bra-side $Y_k$
and the $m{=}0$ bra by the $(1+M)$-dressed versions:
$$
\hat Y_k = R_k X^\dagger(1{+}M)\eta = Y_k + \text{mass\_coeff}\,W_k,\quad
W_k\equiv R_k X^\dagger(M_\text{mass}\eta)\ \ (\textbf{extra multishift solve, RHS }X^\dagger M_\text{mass}\eta),
$$
$$
\hat\Psi=(1{+}M)\eta=\eta+\text{mass\_coeff}\,(M_\text{mass}\eta),\qquad
F=-\tfrac12\mathrm{Re}\Big\{\hat\Psi^\dagger X'(\textstyle\sum Z_k)-\sum_k\big[(X\hat Y_k)^\dagger X'Z_k+\tilde Z_k^\dagger X'\hat Y_k\big]\Big\}.
$$
$\tilde Z_k=XZ_k$ and $Z_k$ are KET-side, untouched. **For massless ($M{=}0$): $\hat Y_k=Y_k$,
$\hat\Psi=\eta$ — identical.** For scalar $m$ it reduces to the old $(1+m^*)$ fold (scalar factors
out of $R_kX^\dagger$). Implementation: in `precalc_grad`, after the $Z_m,Y_m$ solves, when
$\text{mass}\neq0$ solve $W_m=R_m X^\dagger(M_\text{mass}\eta)$ and overwrite `d_Ys[m]`$\to\hat Y_m$
(`d_Ys[m] += mass_coeff*W_m`); store `d_eta_bra`$=\hat\Psi$. Then `grad`/`l1`/`l2`/`l4` use
`d_Ys[m]`($=\hat Y_m$) + `d_eta_bra` and **drop the `conj(mass_coeff)*M_mass` fold entirely**
(the X-precompute `d_XYpre`/`d_XYg` automatically become $X\hat Y_m$ since they are built from the
overwritten `d_Ys[m]`). Validate: L=2 force-vs-FD `~1e-4 -> ~1e-8`, dH-scaling ratio `~4`.

## 5. Fix (SUPERSEDED by 5' -- premult-on-tilde is a no-op for real m)

Implement the massive force EXACTLY as the generalized D.7: form three NEW bra-side copies by a
cheap diagonal $(1+M)$ multiply on the already-solved, X-applied vectors (NO new solve)
$$
\Psi_M \equiv (1{+}M)\,\Psi=(1{+}M)\,\texttt{d\_eta},\quad
\tilde Y_k^M\equiv(1{+}M)\,X Y_m=(1{+}M)\,\texttt{d\_XYs[m]},\quad
\tilde Z_k^M\equiv(1{+}M)\,X Z_m=(1{+}M)\,\texttt{d\_XZs[m]},
$$
and use them as the **bra** in the three force terms:
$$
F^{(m)}=-\tfrac12\mathrm{Re}\Big\{\Psi_M^\dagger X'(\textstyle\sum_k Z_k)
-\sum_k\big[(\tilde Y_k^M)^\dagger X' Z_k+(\tilde Z_k^M)^\dagger X' Y_k\big]\Big\}.
$$
Then **remove** the current ad-hoc `conj(mass_coeff)*M_mass`-on-the-finished-dot fold.

**Why the current fold is wrong (the ordering).** Current term b folds $M_\text{mass}$ onto the
KET `coo Z_m`; by $M_\text{mass}$ Hermiticity that equals $M$ on the bra $\tilde Y_k=XY_m$ — OK for
term b. Current term a folds $M_\text{mass}$ onto `X Z_m` = the bra $\tilde Z_k$ — also looks OK
per-term. The discrepancy is in the **combination** (the $X'=i\,W^{wz}$ factor + the cross-term/sign
structure spoils the naive per-dot $\mathrm{Re}$-equivalence for non-uniform $V$). The robust fix is
NOT to patch slots but to build $\Psi_M,\tilde Y_k^M,\tilde Z_k^M$ explicitly and contract per the
formula above — i.e. premultiply $(1+M)$ on the bra **before** forming each inner product, exactly as
in the scalar case (where $(1+m^*)$ trivially factored out). **Open:** strictly, $(1+M)\tilde Y_k
=(1+M)XR_kX^\dagger\Psi \neq XR_kX^\dagger(1+M)\Psi$ (the exact bra-correction), so this premult-on-
tilde form is exact only if our $C_o$/normalization makes it so (NM asserts it does for this additive
diagonal setting). CONFIRM with NM, then validate by the L=2 force-vs-FD collapsing `~1e-4 -> ~1e-8`.

## 6. Fix + validation plan

1. Confirm the exact $M_\text{mass}$ slot (Sec 5) with NM / re-derivation.
2. Implement in `grad_deviceAsyncLaunch` and `grad_l1`/`grad_l2`/`grad_l4`; **remove** the post-hoc
   `(1+M^*)` = `conj(mass_coeff)*M_mass` fold on the finished inner products. Keep the (correct) $m=0$
   explicit-bra handling. `precalc_grad` likely unchanged (no new RHS) if Sec 5 holds.
3. Validate: `test_diag_mass_l1_claude.cu` L=2 force-vs-FD must drop `~1e-4` $\to$ `~1e-8` (all links,
   real + imag $m$); then `test_hmc_diag_mass_claude.cu` Test 2 dH-scaling must give ratio `~4` at the
   normal `tmax` (the `tmax_sc=0.2` workaround can then be relaxed).
4. Re-run the full 10-phase `tmp_claude.sh` (operator suite) — should stay green; the FD rows are the
   ones that change.

## 7. Code pointers

- Force fold (WRONG slot), default grad pole loop: `includes/overlap_wmass_claude.h` ~`:735-760`
  (`Mp.on_gpuAsync(d_Ms[m], d_XZs[m])` = `M_mass (X Z_m)`; `tmp2reduce[m] -= (i0 + conj(mass_coeff) iM)`).
- `grad_l1` ~`:805-826`; `grad_l2` ~`:874-893`; `grad_l4` ~`:939-958` (same fold, block form via
  `mult_coo_block(M_mass)` on `d_XZg`/`d_CZ`).
- $m=0$ term (CORRECT) ~`:771-780`. `precalc_grad` (builds $Z_m$,$Y_m$) ~`:605-654`.
- Validation drivers: `test_diag_mass_l1_claude.cu` (force-vs-FD), `test_hmc_diag_mass_claude.cu`
  (dH-scaling). Conserved-current $K$: `includes/conserved_current_claude.h` (mass-independent).
