# Krylov-basis GEVP from a single correlator (impl plan)

Source of the idea: M. Wagman, "Lanczos, the transfer matrix, and the signal-to-noise problem",
arXiv:2406.20009 (transfer-matrix reinterpretation of a single correlator). Here we keep the
reinterpretation but replace the Lanczos recursion by the standard GEVP
(Luscher-Wolff, Nucl.Phys.B339 (1990) 222; Blossier et al., arXiv:0902.1265).

## Goal

Extract $E_0$ (and $E_1$) from the region of the effective-mass plot where $m_\text{eff}(t)$ is
still an exponential approach to the plateau (currently fitted with `m0 + A exp(-c1 dt)`, see
`effmass_*_expfit_claude.py`). Instead of modelling the approach, diagonalize it away with a
2x2 GEVP built from the single correlator.

## Construction (i, j = 0, 1 only)

Transfer matrix $T = e^{-aH}$, source state $|\psi\rangle$, $C(t) = \langle\psi|T^t|\psi\rangle$.
Krylov states $|\psi,i\rangle = T^i|\psi\rangle$. Then

$$C_{ij}(t) = \langle\psi,i|\,T^t\,|\psi,j\rangle = C(t+i+j), \qquad i,j\in\{0,1\}$$

i.e. the Hankel matrix

$$
C(t) =
\begin{pmatrix} C(t) & C(t+1) \\ C(t+1) & C(t+2) \end{pmatrix}
$$

GEVP: $C(t)\,v_n = \lambda_n(t,t_0)\,C(t_0)\,v_n$, with
$\lambda_n(t,t_0) = e^{-E_n (t-t_0)}\,[1 + O(e^{-\Delta E (t-t_0)})]$ where, for a 2-state GEVP,
the correction is governed by $E_2 - E_n$ rather than $E_1 - E_0$. That is the whole gain: the
first excited state is removed from the ground-state channel at every $t$, so
$m^\text{GEVP}_{\text{eff},0}(t) = \log[\lambda_0(t,t_0)/\lambda_0(t+1,t_0)]$ should plateau from
much smaller $t$ than the naive $m_\text{eff}(t)$.

Two-state check (exact): if $C(t) = A_0 e^{-E_0 t} + A_1 e^{-E_1 t}$ the 2x2 GEVP is exact for
all $t, t_0$ and returns $e^{-E_0(t-t_0)}$, $e^{-E_1(t-t_0)}$ with no $t$ dependence. This is the
unit test.

Periodic time: for the cosh correlator $C(t) = A(e^{-Et} + e^{-E(N_t-t)})$ the Hankel entry
$C(t+i+j)$ is not literally $\langle\psi,i|T^t|\psi,j\rangle$. Primary target is the small-$t$
region ($t+2 \ll N_t/2$, $N_t = 128$), where the backward term is negligible; a two-state
cosh model can be used as a cross-check later if needed.

Choice of $t_0$: scan $t_0 \in \{1, 2, 3\}$; standard advice $t \le 2 t_0$ is not needed for a
2x2 matrix, but report the $t_0$ dependence.

## Noise / conditioning

$C(t_0)$ is 2x2 and nearly rank-1 when $C(t_0+1)/C(t_0) \approx C(t_0+2)/C(t_0+1)$, i.e. when the
correlator is already single-exponential at $t_0$. That is exactly the case where the GEVP is
not needed; in the exponential-approach region the matrix is well conditioned. Use
`scipy.linalg.eigh(C(t), C(t0))` per jackknife sample; if `eigh` fails (non-PD $C(t_0)$ on a
sample) fall back to `scipy.linalg.eig` and sort by $|\lambda|$.

Errors: config-jackknife of the full pipeline (average $\to$ Hankel $\to$ GEVP $\to$
$m_\text{eff}$), same delete-1 scheme as `effmass_axial_tp_expfit_claude.py:90` (`effmass_jk`).

## Files

- `krylov_gevp_claude.py` (new, `src/production/`): library + driver.
  - `hankel2(C)`: returns `C2[t] = [[C[t], C[t+1]], [C[t+1], C[t+2]]]` for `t = 0 .. Nt-3`.
  - `gevp2(C2, t0)`: returns `lam[n, t]` sorted, `n=0,1`.
  - `meff_gevp(lam)`: `log(lam[:, t] / lam[:, t+1])`.
  - jackknife wrapper reusing `load`, `conn_files`, `C_axial` from
    `effmass_axial_tp_expfit_claude.py` (import or copy; see Q1).
  - synthetic two-state unit test (`--test`): exact recovery of $E_0$, $E_1$.
- `krylov_gevp_<channel>_claude.py` (new): per-channel driver producing PNGs in
  `src/production/figs/` that overlay the naive `meff_acosh` points with the GEVP
  $m_{\text{eff},0}(t)$ and $m_{\text{eff},1}(t)$ for $t_0 \in \{1,2,3\}$.
- `krylov_gevp_masses_claude.md` (new): table of plateau masses vs the existing
  `*_expfit_masses_claude.md` values for the same (Nf, L, g).

No existing file is modified.

## Implementation chunks

1. Core + unit test. Files: `krylov_gevp_claude.py`.
   Hankel, GEVP, effmass, jackknife; synthetic two-state and three-state test
   (three-state: check that the residual $t$ dependence is controlled by $E_2 - E_0$).
2. First real channel: axial tp $\ell=1$ (conn, m-averaged), all (Nf, L, g), since its
   const+exp windows and masses already exist for comparison. Files:
   `krylov_gevp_axial_tp_claude.py`, `figs/krylov_gevp_axial_tp_l1_*.png`.
   Deliverable: side-by-side plateau-mass table vs `axial_tp_l1_expfit_masses_claude.md`.
3. Plateau fits of the GEVP effmass (constant fit over a window chosen from the GEVP curve),
   with jackknife error; write `krylov_gevp_masses_claude.md`. Files: same as chunk 2.
4. Extend to the other channels (scalar PS/FS, vector, glue single-shape correlators) by
   switching the loader. Files: one thin driver per channel.

## Open questions

- Q1. Reuse by `import` from `effmass_axial_tp_expfit_claude.py` (it runs a script body on
  import) or copy the 4 loader functions into `krylov_gevp_claude.py`? Default: copy.
- Q2. Which channel first? Default: axial tp $\ell=1$ (chunk 2) because its reference masses
  exist for all (Nf, L, g).
- Q3. Report $E_1$ from $\lambda_1$ as a physics number, or only use it as a diagnostic?
  Default: plot it, do not tabulate it until it shows a plateau.
