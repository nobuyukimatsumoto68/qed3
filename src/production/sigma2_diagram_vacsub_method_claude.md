# Vacuum subtraction of $\langle\sigma^2\sigma^2\rangle$ (and cross-correlators) by per-diagram analysis

Thorough method note for building a **vacuum-free connected correlator matrix** in the interacting
$0^{++}$ fermionic sector, by (i) decomposing each two-point function into its Wick diagrams, (ii)
identifying which diagrams carry a physical state vs a vacuum constant, and (iii) retaining only the
physical (connected) pieces. This is the input to the GEVP; the variational enhancement (block-Hankel +
rebase) is documented separately in `hankel_rebase_method_claude.md`. Distillation on $S^2\times\mathbb R$,
overlap fermions; driver correlators in `diag_effmass_claude.py::diags_pair` and the `diag_*_claude.py` tools.

## 0. Why this is needed

In the interacting theory the scalar density $\sigma=\bar\psi\Phi\psi$ ($\Phi=$ Y00 $\ell=0$ distillation
smear) has a large condensate. $\langle\sigma^2(t)\sigma^2(0)\rangle$ is then **dominated by a
$t$-independent vacuum constant**, and a naive identity-based GEVP produces a spurious level whose
effmass drifts to $0$ (the condensate/disconnected fluctuation). Subtracting the constant is subtle:
- The **GW contact** $\tfrac12$ (from $D^{-1}+D^{-\dagger}=1$, exact for overlap, interacting included) is
  subtracted on every equal-time perambulator diagonal, `tt[a,a]=tau[a,a]-0.5*I` (`CONTACT=0.5`). This
  fixes the *mean* of $\sigma$ but NOT the disconnected *variance* of $\sigma^2$.
- A single identity subtraction under-removes: the disconnected diagram averages to a **4th moment**
  $\langle(\text{loop})^2(\text{loop})^2\rangle$ while the identity only subtracts $\langle\sigma^2\rangle^2$
  (2nd moment squared) — the mismatch is the leftover constant mode.

**Resolution adopted: vacuum-subtract by DIAGRAM SELECTION, not by fitting/subtracting a constant.**

## 1. The 10 diagrams of $\langle\sigma^2\sigma^2\rangle$ (PS channel)

`diags_pair(Phi, leg, s, t)` returns 10 diagrams (translation-averaged over source $s$, sink $t$; the
full correlator is $2\sum_i W_{10}[i]\,\text{diag}_i$, $W_{10}=[4,2,4,4,2,1,4,1,1,1]$). Building blocks
(source and sink each have two $\sigma$'s; $D_S=\mathrm{Tr}[\Phi\,\tilde\tau]$ is the one-$\sigma$ tadpole,
$D'=\mathrm{Tr}[\Phi\tilde\tau\Phi\tilde\tau]$ the $\sigma^2$ self-loop, $C_S=\mathrm{Tr}[\Phi\tau\Phi\tau]$
the single-$\sigma$ meson propagator, $T_S,V_S,S_S$ the connected 3/4-point loops):

| idx | label | expression | type |
|---|---|---|---|
| 0 | A | $-S_S$ | connected (1 loop, all 4 $\sigma$) |
| 1 | B | $-T_S$ | connected (crossed) |
| 2 | C | $D_S V_S$ | 1 tadpole $\times$ connected 3pt |
| 3 | D | $D_S V_S$ | (= C) |
| 4 | E | $C_S^2$ | **two-meson** ($\langle\sigma\sigma\rangle^2$, 2 loops) |
| 5 | F | $D'D'$ | disconnected (2 self-loops) |
| 6 | G | $-D_S^2 C_S$ | 2 tadpoles $\times$ prop |
| 7 | H | $-D_S^2 D'$ | disconnected |
| 8 | I | $-D_S^2 D'$ | disconnected |
| 9 | J | $D_S^4$ | fully disconnected (glueball-mixing piece; magnitude tiny, ignore) |

**Diagrams carrying a $D_S$ or $D'$ *product-of-self-loops* factor are the vacuum pieces.**

## 2. Per-diagram correlator diagnosis (the key measurement)

Tool: `diag_corr_linear_claude.py` (per-diagram correlator, LINEAR scale, one panel each, plateau/t-sum
subtracted, config jackknife). Companion effmass: `diag_effmass_claude.py`. Findings (Nf2 gsq1.0 L1, 400 cfg):

- **Decay (physical): A, B, E** (and weakly C, D, G). E $=C_S^2$ dominates at small $t$ and decays as the
  two-meson.
- **Flat constants (vacuum): F, H, I** sit at $\sim1.5\times10^{-7}$, do NOT decay, and **dominate the
  TOTAL tail** past $t\sim8$ — this is the constant that drives the identity-GEVP spurious mode. J is tiny.
- The plateau (t-sum) subtraction removes the *bulk* constant, but **F, H, I have a slow residual slope**
  (F: $1.77\to1.47\times10^{-7}$), so a small $\sim0.1$-effmass tail leaks in.

### 2a. Which physical states live where (`diag_BE_vs_2mPS_claude.py`, `diag_sum_vs_mPS_claude.py`)
- **B+E $=C_S^2-T_S$** (direct + crossed, equal weight $W_{10}=2$) = the **two-meson**: effmass plateaus
  $\approx0.70$, i.e. $2m_{PS}+\sim0.05$ — a **weak (repulsive) scattering shift** above threshold
  ($m_{PS}\approx0.32$, $2m_{PS}\approx0.64$; E$=C_S^2$ is the *compact* two-$\sigma$, so a small offset is expected).
- **A+C+D+G** (all $W_{10}=4$) = a **one-meson-rich tower**: effmass $\approx0.46$, between $m_{PS}$ and
  $2m_{PS}$ (an excited-$\sigma$-like scale; ground-$\sigma$ overlap weak in these).
- **A+B+C+D+E+G** = the full **connected** $\sigma^2$: decays cleanly (NO constant), ground $\approx0.48$
  dominated by the light one-meson-rich content, two-meson as a heavier excited piece.

### 2b. Exclusion beats subtraction (`diag_vacsub_sigma2_claude.py`)
Doing "per-diagram plateau subtraction, then a total plateau subtraction, per jackknife" on ALL 10
diagrams does NOT give a clean estimate: F/H/I's slow residual + config variance survive, the effmass
drags to $\sim0.1$, and the correlator hits the noise floor ($C\sim$ its error) by $t\sim12$. **Conclusion:
drop F,H,I,J entirely.** The vacuum-free $\langle\sigma^2\sigma^2\rangle_{\rm conn}=$ A+B+C+D+E+G
(`CONN_IDX=[0,1,2,3,4,6]`).

## 3. Cross-correlator $\langle O_A\,\sigma^2\rangle$ (`diag_OA_sig2_linear_claude.py`)

$O_A=\bar\psi\Phi\tilde\tau\psi$ (bilinear). $\langle O_A(t)\sigma^2(s)\rangle$ has 6 fermion fields $\to$
3 propagators $\to$ 6 Wick pairings in 4 diagram types (multiplicities), sign $=(-1)^{\#\text{loops}}$:

| diagram | mult | expression | role |
|---|---|---|---|
| **Tri** | 2 | $-\mathrm{Tr}[P_A(t)\,\tau(t,s)\Phi(s)\tilde\tau(s,s)\Phi(s)\tau(s,t)]$ | **connected triangle = the physical coupling** |
| Semi | 2 | $+\mathrm{Tr}[P_A\tau(t,s)\Phi\tau(s,t)]\cdot\mathrm{Tr}[\Phi\tilde\tau]$ | $\propto\langle\sigma\rangle\approx0$ (post-contact) |
| OA_Dp | 1 | $+\mathrm{Tr}[P_A\tilde\tau]\cdot\mathrm{Tr}[\Phi\tilde\tau\Phi\tilde\tau]$ | vacuum ($\langle O_A\rangle\langle\sigma^2\rangle$), flat const |
| Disc | 1 | $-\mathrm{Tr}[P_A\tilde\tau]\cdot\mathrm{Tr}[\Phi\tilde\tau]^2$ | $\propto\langle\sigma\rangle^2\approx0$ |

$P_A=\Phi\,\tilde\tau$. **Validation:** $2\,$Tri $=C_{2A}$, the connected cross already used in the GEVP
(matches `gevp_OA_sig22`). So the vacuum-free $\langle O_A\sigma^2\rangle=$ **Tri only** (drop OA_Dp;
Semi/Disc already $\approx0$ because the contact-subtracted $\langle\sigma\rangle$ vanishes).

## 4. The other operators / matrix elements

- Single-$\sigma$ meson propagator $M(t,s)=-\mathrm{Tr}[\Phi(t)\tau(t,s)\Phi(s)\tau(s,t)]$ (= $C_S$),
  effmass $=m_{PS}$; used to build the two-$\sigma$.
- **Time-split two-$\sigma$** $O_{2\sigma}(t)=\sigma_{00}(t)\sigma_{00}(t{+}\delta)$ (`SPLIT`, default 1):
  the genuine four-fermion two-meson interpolator. Its connected pieces are products of $M$:
  $\langle O_{2\sigma}O_{2\sigma}\rangle_c=M(t,s)M(t{+}\delta,s{+}\delta)+M(t,s{+}\delta)M(t{+}\delta,s)$;
  $\langle\sigma^2 O_{2\sigma}\rangle_c=2\,M(t,s)M(t{+}\delta,s)$ — automatically vacuum-free ($M\to0$).
  (A bilinear cannot reach the 4-fermion two-meson — free-theory particle number — hence the split.)
- $O_{22}=\bar\psi Q_2\tilde\tau\psi$ **DROPPED interacting**: gauge broadens the $S^2$ shells past their
  spacing (no 4/8/12 clustering at gsq 0.5/1.0/1.5, checked), so the $\lambda=2$ projector $Q_2$ is
  ill-defined. (Works in the free theory.)

## 5. Assembling the vacuum-free connected GEVP (no identity)

Because every entry is a connected correlator that decays, **no identity operator is needed** (nothing to
subtract), and there is **no vacuum mode**. GEVP energies are invariant under per-operator rescaling, so
the overall diagram normalizations (the $2\sum W_{10}$ convention etc.) do not affect the spectrum.

$$
2\text{-op } \{\sigma^2,O_{2\sigma}\}:\quad
C(t)=\begin{pmatrix} 2\!\sum_{\text{A..G}}\!W_{10}\text{diag} & C_{2s}\\ C_{2s} & C_{ss}\end{pmatrix},\qquad
3\text{-op adds } O_A:\ C_{2A}=2\text{Tri},\ C_{AA}=\text{loop},\ \langle O_A O_{2\sigma}\rangle\approx0.
$$

Config-averaged, jackknife errors; **quote at binsize $\approx$10–16** (autocorrelation, see hankel note).

**Result (Nf2 gsq1.0 L1, 400 cfg, bin 10):** $m_0\approx0.46(1)$ (one-meson-rich ground), $m_1\approx0.644(3)
\approx2m_{PS}$ (two-meson). No vacuum mode; full rank.

## 6. Recipe (summary)

1. Compute all 10 $\langle\sigma^2\sigma^2\rangle$ diagrams; **plot per-diagram linear, plateau-subtracted**
   (`diag_corr_linear_claude.py`) -> confirm A,B,(C,D),E,(G) decay; F,H,I flat; J tiny.
2. Cross-check the physical assignment: B+E $\to2m_{PS}$; A+C+D+G $\to$ one-meson tower
   (`diag_BE_vs_2mPS`, `diag_sum_vs_mPS`).
3. **Keep A,B,C,D,E,G; drop F,H,I,J** (exclusion, not subtraction — `diag_vacsub_sigma2` shows subtraction fails).
4. $\langle O_A\sigma^2\rangle$ -> **triangle only** ($C_{2A}$); $O_{2\sigma}$ from $M$ products.
5. Assemble the connected matrix, **no identity**; GEVP; binsize-corrected errors.

## 7. Files

- Diagram diagnosis: `diag_effmass_claude.py` (10-diagram effmass + `diags_pair`),
  `diag_corr_linear_claude.py` (per-diagram linear corr, panels), `diag_BE_vs_2mPS_claude.py`,
  `diag_sum_vs_mPS_claude.py` (`DIAGS=` any subset), `diag_vacsub_sigma2_claude.py` (double-subtraction test),
  `diag_OA_sig2_linear_claude.py` (⟨$O_A\sigma^2$⟩ diagrams).
- GEVP driver: `gevp_twosigma_interacting_claude.py` (`CONN=1`, `NOOA` for 2-op).
- Companion: `hankel_rebase_method_claude.md` (variational enhancement).

Refs: distillation Peardon et al. 0905.2160; mixing Chester–Pufu 1603.05582.
