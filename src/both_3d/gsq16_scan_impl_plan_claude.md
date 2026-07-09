# Add gsq=16.0 to the L=1 massless gsq scan -- plan

## Goal
Extend the massless L=1 coupling scan (gsq = 1, 2, 4, 12 done/running) with **gsq=16.0**,
Nf=2/4/6, to strengthen the strong-coupling end of the scan.

## What's identical to the rest of the scan (no code changes)
- Source: `hmc_L1_claude.cu` (unchanged) -- npole=21, window k=0.001 (delta~1.5e-5),
  nsteps=10 all Nf, tmax=1.9, at=0.2, Nt=128, massless, k_ckpoint=1, k_ckpoint_rng=10.
- gsq is a runtime arg -> just `./hmc_L1_claude.o 16.0 <Nf> 1.0`. Fresh dirs
  `Nf{2,4,6}_gsq16.000000at0.200000nu01.000000nt128L1/` (no collision, verified).
- Only new file: `run_L1_gsq16_claude.sh` (copy of run_L1_gsq12 with GSQ=16.0 + the
  MPS-auto-start block already standard).

## Cost / risk (the real content)
- gsq=16 is STRONGER coupling than 12 -> rougher gauge fields -> more inner-CG iterations
  (TOL_INNER=1e-9) -> SLOWER per traj than gsq12. Reference (gsq12, but measured while MPS
  was down/degraded): Nf2 ~647s, Nf4 ~711s, Nf6 ~1113s. Post-reboot with MPS up these are
  ~5-6x faster, so a healthy gsq16 3-pack is plausibly ~150-400 s/traj/stream -> to kmax=320
  from k=0, order ~1-3 days for the pack (Nf6 the long pole).
- FREEZE RISK is higher at strong coupling (more near-zero modes). npole21/window0.001 is the
  robust fix and held for gsq12, but watch via `/check-configs`; recovery = rollback + finer
  nsteps (per skill playbook).
- MPS: must be UP (run script auto-starts it). 3/GPU vs 2/GPU: strong coupling is CG-heavy so
  MPS payoff is modest; 2-3 per GPU is fine, judge by ckpoint-mtime throughput after launch.

## Open questions (resolve before launching)
1. **kmax**: 320 (match the scan) or a shorter cap (e.g. 160) given the strong-coupling cost?
2. **When / which GPU**: GPU1 is busy with gsq=4 (finishes in ~8h); GPU0 has the ylm-scalar
   job (other session). Options: (a) queue gsq16 on GPU1 after gsq=4 completes; (b) start now
   packed with gsq=4 on GPU1 (more contention); (c) wait for GPU0.
3. Packing: 3/GPU (Nf2/4/6 together) or the {Nf2->Nf4}+{Nf6} 2-slot grouping used for gsq=4?
4. Is gsq=16 the last point, or are more couplings coming (affects whether to templatize the
   run script to take gsq as an arg)?
