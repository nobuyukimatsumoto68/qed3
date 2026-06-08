# Audit: jj_corr_block_t_claude.cu vs PDF (v2-11) -- operators, GW factors, leg structures

Verification that the connected/disc estimators in `jj_corr_block_t_claude.cu` implement the PDF formulas
(Sec. 3.3.2 massless, 3.3.3 m_F, 3.3.4 m_P).  Focus: D, Dm, Dtil and the $(1-D_\text{ov})$ insertions.

**VERDICT: no bug found.**  The operators, GW factors, and all measured leg structures match the PDF.
One modeling choice (m_F axial = massless-limit legs) to confirm; it does NOT affect the massless run.

## Operators (code lines 255-257; overlap_wmass_claude.h)

| symbol | PDF | code | verdict |
|---|---|---|---|
| $D=D_\text{ov}$ | massless overlap | `Fermion D(DW, 0, 11)` | OK |
| $D_m=D_\text{ov}+m$ | Eq. 3.61 | `Fermion Dm(DW, valence_mass, 11)` | OK |
| $\tilde D_{m_P}=D_\text{ov}+\frac{m_P}{1-m_P}$ | Eq. 3.64 | `Fermion Dtil(DW, valence_mass/(1-valence_mass), 11)` | OK |

`mult` adds `+m v` (line 373/433); `adj` adds `+m^* v` (line 405/456) = $(D_\text{ov}+m)^\dagger=D_\text{ov}^\dagger+m^*$. OK.
Solver convention: $D_\text{ov}+m$ is NORMAL, so `solve_sq` gives $(X^\dagger X)^{-1}X^\dagger=X^{-1}$ and
$(X^\dagger X)^{-1}X=X^{-\dagger}$ -- used consistently below.

## $(1-D_\text{ov})$ GW insertions (op_oneMinusD, lines 274-276)

`op_oneMinusD = 1 - D_ov` (identity + (-D_ov)).  Applied to the **source** as
$\rho=(1-D_\text{ov})K^\dagger(t_0)\eta$ -> after the source dagger it appears as $(1-D_\text{ov}^\dagger)$ in
the trace; and to the **forward leg** $\chi=(1-D_\text{ov})\phi'_A$ -> appears as $(1-D_\text{ov})$ at the
sink.  Matches the PDF's $(1-D_\text{ov}^\dagger)$ on the $\eta$/$t_0$ side and $(1-D_\text{ov})$ on the
$\xi$/$t$ side (Eqs. 3.55, 3.69).  OK.

## MASSLESS (the immediate Nf=2,4 run) -- Eqs. 3.52-3.56

- **Vector $C_{V++,c}$ (3.52)** $=\mathrm{tr}[K(t_0)D_\text{ov}^{-1}K(t)D_\text{ov}^{-1}]$.
  Code: $\phi'=D_m^{-1}\eta$, source $\psi=D_m^{-\dagger}K^\dagger(t_0)\eta$ (op_Dm + solve_sq), sink
  $K(t)\phi'$; dot $\psi^\dagger(K\phi')=\eta^\dagger K(t_0)D_m^{-1}K(t)D_m^{-1}\eta$.  At $m=0$ ($D_m=D_\text{ov}$). **OK.**
- **Disc $C_{V++,d}$ (3.53)** $=\mathrm{tr}(D_\text{ov}^{-1}K)\,\mathrm{tr}(D_\text{ov}^{-1}K)$: $J(t)=\eta^\dagger K(t)\phi'$. **OK.**
- **$C_{V--}=C_{V++}^*$ (3.54):** written via `write_corr_conj` (non-parity). **OK.**
- **Axial $C_{A+-}$ (3.55)** $=\mathrm{tr}[K(t_0)(1-D_\text{ov}^\dagger)D_\text{ov}^{\dagger-1}\,K^\dagger(t)(1-D_\text{ov})D_\text{ov}^{-1}]$.
  Code: forward $\phi'_A=D_m^{-1}\eta$, $\chi=(1-D_\text{ov})\phi'_A$; source
  $\psi=D_m^{-1}(1-D_\text{ov})K^\dagger(t_0)\eta$ (op_DmH + solve_sq) -> $\psi^\dagger=\eta^\dagger
  K(t_0)(1-D_\text{ov}^\dagger)D_m^{-\dagger}$; sink $K^\dagger(t)\chi$ (apply_k_**dag**).
  Dot -> $\mathrm{tr}[K(t_0)(1-D_\text{ov}^\dagger)D_m^{-\dagger}K^\dagger(t)(1-D_\text{ov})D_m^{-1}]$, $m=0$. **OK.**
  - **Sink dagger note (the one thing that looks odd, but is correct):** the code applies $K^\dagger$ at the
    sink.  $C_{A+-}=\langle J_{A,+}(t_0)J_{A,-}(t)\rangle$ with $J_{A,+}\sim K$, $J_{A,-}\sim K^\dagger$, so
    the **$-$ current at the sink uses $K^\dagger$** -- consistent.  The OCR of (3.55) shows the second
    kernel as un-daggered $K(t)$, but (3.56) PRESERVES the analogous $K^\dagger(t_0)$, confirming the
    $+/-$ kernel convention; (3.55)'s sink $\dagger$ is just dropped in the text rendering.  Code matches.

## m_P (parity, later runs) -- Eqs. 3.65-3.70

- **Vector $C_{V++}$ (3.65)** uses $D_{m_P}=D_\text{ov}+i m_P$: code (++) uses op_Dm with `valence_mass=i m_P`. **OK.**
- **Vector $C_{V--}$ (3.67)** $=\mathrm{tr}[K^\dagger(t_0)\tilde D_{m_P}^{\dagger-1}K^\dagger(t)\tilde D_{m_P}^{\dagger-1}]$:
  code (--) builds $\psi=\tilde D^{-1}K(t_0)\eta$ (op_tilDmH+solve_sq), sink $K^\dagger(t)\,\mathrm{phimm}$,
  $\mathrm{phimm}=\tilde D^{-\dagger}\eta$ -> $\mathrm{tr}[K^\dagger(t_0)\tilde D^{\dagger-1}K^\dagger(t)\tilde D^{\dagger-1}]$. **OK.**
- **Axial $C_{A+-}$ (3.69)** $=\mathrm{tr}[K(t_0)(1-D_\text{ov}^\dagger)\tilde D_{m_P}^{\dagger-1}K^\dagger(t)(1-D_\text{ov})D_{m_P}^{-1}]$:
  code axial parity uses op_tilDmH (X_sink=tilde) for the source solve and $\phi'_A=D_m^{-1}\eta$ -> matches
  ($\tilde D^{\dagger-1}$ on $t_0$, $D_{m_P}^{-1}$ on $t$). **OK.**
- The $(1+m_P)^{-1}$ of (3.63) is cancelled by the adjoint-kernel correction (PDF note below 3.64); code
  uses $\tilde D$ directly (no explicit $(1+m_P)$).  Consistent. **OK.**

## m_F (flavor) -- Eqs. 3.59-3.61 + note

- **Vector:** "(3.52-3.54) with $D_\text{ov}\to D_{m_F}$" -> code (++) uses op_Dm = $D_\text{ov}+m_F$. **OK.**
- **Axial:** PDF gives NO explicit m_F axial formula.  Code uses the **massless** $D_\text{ov}$ legs
  (op_DH/op_Dsq; `flavor` branch) -- the "no conserved axial current at $m_F\neq0$ => massless-limit"
  modeling choice (comment line 251-203).  **CONFIRM this is intended** (not derivable from the PDF as
  written).  Irrelevant to the massless run.

## Conclusion

For the **massless Nf=2,4 run**: operators, GW factors, and the vector + axial massless estimators all
match the PDF -- safe to run.  The m_P legs match 3.65-3.70.  Only open item: confirm the m_F-axial
massless-limit modeling (does not affect massless).
