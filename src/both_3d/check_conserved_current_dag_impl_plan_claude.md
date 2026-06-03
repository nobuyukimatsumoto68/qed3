# Implementation Plan: Adjoint Kernel $K^{wz\dagger}$ and its Check Program

_Created: 2026-06-03_

---

## 1. Physics / Goal

Implement the adjoint action of the exactly conserved vector current kernel,
$(K^{wz})^\dagger$, as a new method `apply_k_dag` in `conserved_current_claude.h`,
and verify it with a check program structured exactly like the existing
`check_conserved_current_claude.cu`.

### 1.1 Structure of the forward kernel $K$

From `conserved_current_claude.h:71-139`, with the operators

- $X = D_W/\lambda_\text{max}$ (apply via `M_DW`, scale $1/\lambda_\text{max}$);
  $X^\dagger = D_W^\dagger/\lambda_\text{max}$ (apply via `M_DWH`).
- $R_m = \left(X X^\dagger - (k^2/c_m)\,I\right)^{-1}$, built as
  `(1/lmax^2) M_DW M_DWH - k^2/cp[m]`. Since $X X^\dagger$ is Hermitian and the
  shift is real, $R_m^\dagger = R_m$.
- $C \equiv$ `D.C` (real), $A_m \equiv$ `D.A[m]` (real).
- $W \equiv W^{wz}$ (built by `build_W`), $W^\dagger$ (built inline in `apply_k`).
- $S \equiv I + \sum_{m\ge 1} A_m R_m$ (Hermitian, $S^\dagger = S$).

the forward kernel is

$$
K = C\,W\,S \;-\; C\sum_{m\ge 1} A_m\, X R_m \left(X^\dagger W - W^\dagger X\right) R_m .
$$

Here $Z_m = R_m\,\xi$, $Z_0 = S\,\xi$, $u_m = (X^\dagger W - W^\dagger X) R_m\,\xi$,
Term A $= C\,W Z_0$, Term B $= -C\sum_m A_m X R_m u_m$.

### 1.2 The adjoint $K^\dagger$

Using $C, A_m\in\mathbb{R}$, $R_m^\dagger = R_m$, $S^\dagger = S$, and
$(X^\dagger W - W^\dagger X)^\dagger = W^\dagger X - X^\dagger W$:

$$
K^\dagger = C\,S\,W^\dagger \;-\; C\sum_{m\ge 1} A_m\, R_m \left(W^\dagger X - X^\dagger W\right) R_m\, X^\dagger .
$$

This is the same shape as $K$ with $W \leftrightarrow W^\dagger$ swapped, $X
\leftrightarrow X^\dagger$ swapped, the operator order reversed, and the inner
antisymmetric block sign-flipped. Implementation:

**Term A$^\dagger$** $= C\,S\,W^\dagger\,\xi$:
1. $w_0 = W^\dagger\xi$ (`coo_WH(w0, xi)`).
2. $Z_0' = S\,w_0 = w_0 + \sum_m A_m R_m w_0$ (one $R_m$ solve per pole on $w_0$).
3. result $= C\,Z_0'$.

**Term B$^\dagger$** $= -C\sum_m A_m\, R_m (W^\dagger X - X^\dagger W) R_m X^\dagger\,\xi$:
1. $p = X^\dagger\xi = (1/\lambda_\text{max})\,D_W^\dagger\,\xi$.
2. For each pole $m$:
   - $q_m = R_m\,p$ (solve).
   - $a = W^\dagger (X q_m)$, with $X q_m = (1/\lambda_\text{max})\,D_W q_m$, then `coo_WH`.
   - $b = X^\dagger (W q_m)$, with $W q_m =$ `coo_W(q_m)`, then $(1/\lambda_\text{max})\,D_W^\dagger$.
   - $r_m = a - b = (W^\dagger X - X^\dagger W) q_m$.
   - $s_m = R_m\,r_m$ (solve).
   - result $\mathrel{-}= C A_m\, s_m$.

Cost: Term A is one $R_m$ solve per pole on $w_0$; Term B is two $R_m$ solves per
pole ($q_m$ then $s_m$) -- same order as the forward `apply_k`.

### 1.3 Fundamental correctness relation (primary check)

For every link $\ell$ and random $\eta,\phi$:

$$
\langle \eta,\; K^\ell\,\phi\rangle \;=\; \langle K^{\ell\dagger}\eta,\; \phi\rangle .
$$

This validates `apply_k_dag` directly against the already-trusted `apply_k`, with
no new physics assumptions. It is the analog of, and replaces, the force / Theta
checks for the forward kernel.

### 1.4 Mirrored identity checks (analogs of the $K$ battery)

The forward check program verifies (in order): force cross-check, $[X,\Theta]$
check, $[D_\text{ov},\Theta]$ check, $\mathrm{tr}(K)$, and the Ward identity
(Chunks 5b/6). The adjoint analogs are:

- **$\mathrm{tr}(K^\dagger)$:** $\mathrm{tr}(K^{\ell\dagger}) = \overline{\mathrm{tr}(K^\ell)} = 0$.
- **$D_\text{ov}^\dagger$ commutator:** since $\sum_z K^{wz} = [D_\text{ov},\Theta_w]$,
  taking the adjoint gives
  $$
  \sum_z K^{wz\dagger} = [D_\text{ov},\Theta_w]^\dagger = -[D_\text{ov}^\dagger,\Theta_w],
  $$
  i.e. $[D_\text{ov}^\dagger,\Theta]\,\xi = -\sum_\ell \delta\theta_\ell\, K^{\ell\dagger}\xi$.
- **Ward identity:** $\sum_{z\sim w}\mathrm{tr}(K^{wz\dagger} D_\text{ov}^{-1}) =
  \overline{\sum_{z\sim w}\mathrm{tr}(D_\text{ov}^{-\dagger} K^{wz})}$; both Re and Im
  sum to zero config-by-config (Ginsparg-Wilson). Estimator mirrors Chunk 5b/6.

---

## 2. Files

### Modified
- `src/both_3d/includes/conserved_current_claude.h` -- add `apply_k_dag` and
  spatial wrapper `apply_k_dag_wz` next to `apply_k` / `apply_k_wz`.

### Created
- `src/both_3d/check_conserved_current_dag_claude.cu` -- standalone check program,
  copy of `check_conserved_current_claude.cu` with `apply_k` -> `apply_k_dag` in
  the identity checks, plus the new adjoint-consistency block (1.3).

### Read-only references
- `src/both_3d/check_conserved_current_claude.cu` -- structure to mirror.
- `src/both_3d/includes/conserved_current_claude.h` -- forward kernel.

---

## 3. Ordered Implementation Chunks

### Chunk D1 -- `apply_k_dag` in the header

`apply_k_dag` is a **reordering of the same primitives** used by `apply_k`
(`conserved_current_claude.h:71-139`); no new kernels, operator types, COO
formats, sign tables, or device buffers. The subtle hand-coded pieces are reused
verbatim:

- $W$/$W^\dagger$ COO construction including the manual adjoint
  (`-Im/lM, -Re/lM` swap + index transpose, lines 96-105) -- copied unchanged.
- $R_m$ solve operator `(1/lmax^2) M_DW M_DWH - k^2/cp[m]`, the $X$ and
  $X^\dagger$ MatPolys, scratch `d_Zs`, `d_tmp1/2/3` -- all reused.

**Term-by-term mapping** of $K^\dagger = C\,S\,W^\dagger -
C\sum_m A_m R_m(W^\dagger X - X^\dagger W)\,R_m X^\dagger$:

| Piece | vs `apply_k` | New code |
|---|---|---|
| Step 1: seed $Y_m = R_m X^\dagger\xi$ | Step 1 (73-90), rhs $X^\dagger\xi$ instead of $\xi$ (one extra `M_DWH` apply + $1/\lambda_\text{max}$); drop the $Z_0$ accumulation line | ~1 line changed |
| Step 2: build $W$, $W^\dagger$ | identical (92-105) | copy verbatim |
| Term A$^\dagger$ $= C\,S\,W^\dagger\xi$ | new shape: `coo_WH(w0, xi)`; result $= C\,w_0$; per pole $m$: solve $R_m w_0 \to t$, result $\mathrel{+}= C A_m t$ | ~8 lines, existing primitives |
| Term B$^\dagger$ | same 5 primitives as Term B (117-137), reordered: $X^\dagger$ now seeds ($Y_m$), no trailing $X$ apply; per pole: $a = W^\dagger(X Y_m)$, $b = X^\dagger(W Y_m)$, $r = a-b$, solve $R_m r \to s$, result $\mathrel{-}= C A_m s$ | ~18 lines, existing primitives |

Total ~45 lines vs `apply_k`'s ~70; the $W^\dagger$ COO (most sign-bug-prone
piece) is reused identically. The Step-1 seed $Y_m = R_m X^\dagger\xi$ is the
same solve the force precalc (`precalc_grad_deviceAsyncLaunch`, overlap.h:414,
`d_Ys`) already uses, giving independent corroboration. Add `apply_k_dag_wz`
spatial wrapper mirroring `apply_k_wz`.

**Note (sign):** $W^\dagger X - X^\dagger W = -(X^\dagger W - W^\dagger X)$, the
negative of the Term B antisymmetric block in `apply_k`; verify this sign when
forming $r = a-b$.

**Files:** `conserved_current_claude.h`

### Chunk D2 -- Adjoint-consistency check (primary)
New check program boilerplate (copy Chunks 1-3 of the forward file). Add the block:
for random $\eta,\phi$ and every spatial + temporal link, compute
$\langle\eta, K^\ell\phi\rangle$ via `apply_k` and $\langle K^{\ell\dagger}\eta,\phi\rangle$
via `apply_k_dag`; report $\max_\ell$ absolute difference; assert $< 10^{-7}$
(free-field is exact $\sim 10^{-15}$; Gaussian floors at the inner-solver
tolerance `TOL_INNER`$=10^{-9}$ compounded over the $R_m$ solves).
**Files:** `check_conserved_current_dag_claude.cu`

### Chunk D3 -- $[D_\text{ov}^\dagger,\Theta]$ commutator check
Adjoint analog of the forward `[D_ov,Theta]` check, using the relation in 1.4:
$$[D_\text{ov}^\dagger,\Theta]\,\xi = -\sum_\ell \delta\theta_\ell\, K^{\ell\dagger}\xi.$$
With $D_\text{ov}^\dagger$ applied via `adj_deviceAsyncLaunch` (wrapped as `op_DH`):
- $\chi_1 = D_\text{ov}^\dagger(\Theta\xi) - \Theta(D_\text{ov}^\dagger\xi)$.
- $\chi_2 = -\sum_\ell \delta\theta_\ell\, K^{\ell\dagger}\xi$ (spatial:
  $\delta\theta = \theta_w - \theta_z$; temporal: $\theta_{(s,ix)} - \theta_{(s+1,ix)}$;
  **note the overall minus sign**).
Include the $\Theta = c\,I$ sanity case ($\chi_1=\chi_2=0$) and a random-$\Theta$
case; verify $\Vert\chi_1-\chi_2\Vert_\infty$ is small (assert $< 10^{-6}$, set by
the overlap pole solves inside $D_\text{ov}^\dagger$ and `apply_k_dag`).
**Files:** `check_conserved_current_dag_claude.cu`

### Chunk D4 -- $\mathrm{tr}(K^\dagger)$ check
Mirror the `run_trK` block with `apply_k_dag`; exact tr via basis loop and
stochastic $\eta^\dagger K^\dagger \eta$; both $\to 0$.
**Files:** `check_conserved_current_dag_claude.cu`

### Chunk D5 -- Ward identity (Chunks 5b/6 analogs)
Mirror Chunk 5b (free-field) and Chunk 6 (Gaussian $U$) with `apply_k_dag`;
exact basis-vector reference + stochastic estimator; report
`Re_analy/Re_stoch/err_Re/sigma_Re` and Im counterparts; divergence $\to 0$.
**Files:** `check_conserved_current_dag_claude.cu`

---

## 4. Resolved Decisions (2026-06-03)

1. **Scope:** adjoint-consistency (D2) + $[D_\text{ov}^\dagger,\Theta]$ commutator
   (D3) + $\mathrm{tr}(K^\dagger)$ (D4) + Ward identity (D5). (Commutator check
   re-added 2026-06-03 at user request.)
2. **Layout:** new separate file `check_conserved_current_dag_claude.cu`.
3. **Method:** new `apply_k_dag` (+ `apply_k_dag_wz`) in `conserved_current_claude.h`.
