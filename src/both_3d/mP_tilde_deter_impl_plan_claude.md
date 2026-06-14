# Impl plan -- tilde-D support for the m_P deterministic ground truth

Goal: make the parity-asymmetric ($\tilde D_{m_P}$) channels of the FREE deterministic ground truth
computable at $m_P$, so the diluted stochastic $m_P$ run can be validated for vector `Vmm` and axial `Apm`
(currently only vector `Vpp` + condensate are checkable). Ref: qed3int_v2-14 Eqs. 3.60-3.67.

## What's already there
- `jj_propagator_deter_claude.cu` ALREADY saves `Dtil_inv` $=\tilde D_{m_P}^{-1}$ at parity
  (`:281,395,411`), alongside `Dm_inv` $=D_{m_P}^{-1}$ and `Dov`.  $\tilde D_{m_P}=D_{ov}+m_P/(1-m_P)$.
- `condensate_deter_free_claude.cu` already handles $m_P$ fully.

## Key identities (no extra solves)
The contraction programs hold the dense kernel $K$ (and $D_m^{-1}$) as row-major $N\times N$.
- $K^\dagger$ (the "-" current kernel) = CONJUGATE-TRANSPOSE of the dense $K$: $(K^\dagger)_{ij}=\overline{K_{ji}}$.
- $\tilde D_{m_P}^{-\dagger}$ = conjugate-transpose of the loaded `Dtil_inv`.
So everything is a host O(N^2) conj-transpose + the existing `matmul_A`/`conn_shift`/`trace`; NO new solves.

## Vector exact-K (Eqs. 3.63/3.65) -- `jj_exact_diag_deter_free_claude.cu` (single-insertion path)
Vpp (unchanged): $A_{++}=K\cdot D_{m_P}^{-1}$; `Vpp=conn_shift(A_{++})`, `disc=tr(A_{++})`.
Vmm (NEW, parity): $A_{--}=K^\dagger\cdot\tilde D_{m_P}^{-\dagger}$; `Vmm=conn_shift(A_{--})`, `disc_mm=tr(A_{--})`.
- Load `Dtil_inv` at parity; form `Ptil = conjT(Dtil_inv)`.  Per `which` (tp/sp): `Kdag = conjT(K)`,
  `A_mm = matmul_A(Kdag, Ptil)`, write `Vmm = conn_shift(A_mm)` (NOT conj), `disc/<proj>/Jmm`.
- Replace the `if(!parity) write Vmm=conj` with: non-parity -> conj(Vpp); parity -> the A_mm result.

## Axial exact-K (Eq. 3.67) -- `jj_exact_axial_deter_free_claude.cu`
$C_{A+-}=\mathrm{tr}[K(t0)(1-D_{ov}^\dagger)\,\tilde D_{m_P}^{-\dagger}\,K^\dagger(t)(1-D_{ov})\,D_{m_P}^{-1}]$.
The massless/$m_F$ code already builds $B_{src}=K(1-D_{ov}^\dagger)P_{src}$ and $B_{snk}=K^\dagger(1-D_{ov})P_{snk}$
with $P_{src}=P_{snk}=D_m^{-1}$.  For $m_P$ the SINK-leg propagator becomes $\tilde D_{m_P}^{-\dagger}$:
set $P_{snk}=\mathrm{conjT}($`Dtil_inv`$)$ (load it at parity), $P_{src}=D_{m_P}^{-1}$ unchanged.  (Confirm
which leg carries tilde by matching the dilute: dilute parity axial sink uses `op_tilDmH`/`blk_*_Dtil`,
forward uses `op_DmH`/Dm -> sink = tilde.)  Single complex channel `Apm`; no Vmm.

## KEY: (1+m_P)^{-1} factor -- bare currents keep it, the conserved current cancels it
PDF p.18 (after Eq. 3.62): for the CURRENT correlators the adjoint-KERNEL correction CANCELS the
(1+m_P)^{-1} of the propagator Eq. (3.61) -> the exact-K Eqs. (3.65/3.67) have NO (1+m_P)^{-1}.  The BARE
LOCAL current (and the condensate bilinears) have NO such kernel -> the factor is NOT cancelled and the
propagator rule (3.61) applies DIRECTLY, INCLUDING (1+m_P)^{-1}.  So:
- exact-K vector Vmm / axial: tilde D^{-dag}, NO (1+m_P)^{-1}.   [chunks 1,2]
- local vector Vmm: dilute SKIPS it at parity (`:892`) -> nothing to compare (skip).
- local AXIAL + condensate: tilde D^{-dag} WITH (1+m_P)^{-1}.    [chunk 3 + condensate already]

## Local axial -- tilde t0-leg + (1+m_P)^{-1} in BOTH dilute and determ (user-requested fix)
The dilute local-axial t0-leg was D_m (unspecialized for parity).  FIXED both sides:
- `jj_corr_dilute_claude.cu`: at parity the psi_A solve uses `op_tilDmH`+`blk_tp_Dtil` (tilde D^{-1}(sigma eta))
  instead of `op_DmH`+`blk_tp_Dm`; CsA *= (1+m_P)^{-1} at output.
- `jj_local_axial_deter_claude.cu`: `tr_WPWP_axial` generalized to (E0,Et,P0,P) -- t0-leg daggers P0;
  at parity P0=Dtil_inv (-> tilde D^{-dag}), P=Dm_inv, and Cpp/Cyl *= (1+m_P)^{-1}.

## Chunks
1. `jj_exact_diag_deter_free` Vmm (vector exact-K -).            [DONE]
2. `jj_exact_axial_deter_free` tilde t0-leg (axial exact-K).     [DONE]
3. local axial tilde + (1+m_P)^{-1} in dilute + jj_local_axial_deter.  [DONE]
4. `run_determ_mP_claude.sh`: all 6 steps enabled (exact axial + local axial now tilde-correct). [DONE]

## Open questions
- OQ1 (disc at m_P): exact_diag now writes `disc/<proj>/Jmm = tr(A_mm)` (the - disconnected).  DONE.
- OQ2 (tilde leg): CONFIRMED tilde is the t0-leg (B_src/Pdag) in the axial, matching the dilute parity
  block (forward `op_DmH` = D_m, sink `op_tilDmH` = tilde; conn_shift2 puts B_src at t0).
- OQ3 (local axial physics): local axial $m_P$ is $D_m$ in both dilute+determ (consistent, not tilde).
  Make it tilde-correct?  Touches the PRODUCTION dilute local-axial block -> ask the user.
