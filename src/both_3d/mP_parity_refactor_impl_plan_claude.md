# Impl plan -- fix the m_P (parity) blocks in jj_corr_dilute_claude.cu (audit I1/I2/I3/I4)

Audit: jj_corr_dilute_audit_claude.md. The parity `(--)` blocks crash under dilution (write inside the
pattern loop -> duplicate HDF5 key on pattern 2) and are in a different (per-origin) frame than `Vpp`.
Fix: make ALL parity outputs SUPERPOSED + hit-scope-accumulated + written ONCE after the dilution loop,
mirroring the `(++)` path; and ADD the local vector `Vmm` for m_P (I3).

## (1+m_P) power bookkeeping (BARE currents keep it; conserved K cancels it)
Each BACKWARD leg <eta xi^dag> = (1+m_P)^{-1} tilde D^{-dag} (Eq. 3.61); forward <xi eta^dag>=D_m^{-1}
(none).  Conserved-current kernel corrections cancel ALL factors (Eqs. 3.65/3.67 -> none).  Bare currents:
- local V++ : 2 forward -> no factor (current).
- local V-- : 2 BACKWARD -> **(1+m_P)^{-2}**   [NEW, I3].
- local A+- : 1 fwd + 1 bwd -> (1+m_P)^{-1}     [done].
- condensate xidag_eta / xidag_1mDdag_eta : 1 bwd -> (1+m_P)^{-1}  [done].
- exact-K V--/axial : kernel-cancelled -> no factor  [done in deter; dilute V-- below].

## Target output frame (parity), all written ONCE after the dilution loop
- h0/tp/Vmm = Ctp_mm (superposed C(t)=G(t)+G(t-Nt/2));  h0/sp/Vmm = Csp_mm.   (exact-K, NO factor.)
- h0/disc/tp/Jtil = JtpT_til;  h0/disc/sp/Jtil = JspT_til.
- h0/s{c}/Vmm = (1+m_P)^{-2} * Cs_mm[c]  (local vector, NEW).
(Vpp / s{c}/Vpp / axial / condensate already write once after the loop.)

## Estimators (superposed, mirror the (++) path)
- exact V-- tp: rho_mm[n]=sum_b K(n,t0_b)eta (dag=false); psi_mm[n]=tilde D^{-1} rho_mm[n]
  (op_tilDmH + blk_tp_Dtil); sink K^dag(n,t)phimm (phimm=tilde D^{-dag}eta, already solved at :602);
  Ctp_mm[t] += w_tp[n] psi_mm[n].dag(kphi)  (absolute t).  = tr[K^dag(t0) tilde D^{-dag} K^dag(t) tilde D^{-dag}] (Eq.3.65).
- exact V-- sp: analogous over links -> Csp_mm.
- disc Jtil: tilphi=tilde D^{-1}eta (own solve, :729); JtpT_til[t]+=w_tp[n](K(n,t)tilphi)^dag eta; JspT_til over links.
- local V-- (NEW): psi_mm[n]=tilde D^{-1}(sum_b sigma_c(n,t0_b)eta); sink chi_mm=sigma_c phimm;
  Cs_mm[c][t] += w_site[n] psi_mm[n].dag(chi_mm) at (t,n).  -> (1+m_P)^{-2} at output.

## Chunks
1. Add hit-scope accumulators (after :512): Ctp_mm,Csp_mm,JtpT_til,JspT_til (Nt vectors); Cs_mm[3].
2. Refactor `(--)` tp (:599-653): summed source, accumulate Ctp_mm, NO write; KEEP the condensate
   block (:608-610) and phimm solve; drop the inert N_ELL ylm bits; rename local shadow buffer (I4).
3. Refactor `(--)` sp (:687-720): summed source, accumulate Csp_mm, NO write (I4 rename).
4. Refactor disc `(--)` (:726-757): accumulate JtpT_til/JspT_til, NO write (drop inert ylm Jtil).
5. ADD local V-- in the local block (:765-789), parity branch -> Cs_mm.
6. Output (after loop): parity writes h0/{tp,sp}/Vmm, h0/disc/{tp,sp}/Jtil, h0/s{c}/Vmm*(1+m_P)^{-2}.
7. Cosmetic: I5 eq# comments (:290-291 3.60->3.59, 3.63->3.62); I7 one-line note at the flavor flag.

## Notes / risks
- buffers reused: psi_tp[n] (summed, ITP(n,0)), rho, hblk, kphi, kblk, d_sinkvec, phimm, tilphi -- same as
  the (++)/disc paths; phimm valid from :602 through the local block (overwritten only by the conn axial at :833).
- mass-independent K cache + non-parity path UNCHANGED -> massless/m_F production byte-identical.
- After: rebuild jj_corr_dilute_mP.o (run_jj_dilute_valid_mP) -- m_P now runs to completion, Vpp/Vmm same frame.
- I8 (TOL_OUTER) left as-is (separate; flag only).
