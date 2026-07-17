# REDO campaign -- runtime / cost estimate (2026-07-15)

Estimated GPU-hours and wall-clock for the corrected-gsq bug-fixed redo (affine, `qos=opp`, MPS 2 streams/GPU,
4h jobs). See `params_L1L2_claude.md` / `params_massive_claude.md` for the physics scope. Decision (NM):
**apply all 13 pair-jobs on `opp`**.

## Measured anchor
- **L1 Nf2 (Hasenbusch, weak gsq, packed 2/GPU): ~35 s/traj** -- from the smoke jobs (`# HMC : ~35 sec`;
  33 traj in 19.5 min). This is the only DIRECTLY measured redo number.
- All other rows below are EXTRAPOLATED, cross-checked against the OLD strong-coupling ($g^2=12,16$) massless
  logs (L2 Nf2/4/6 = 393 / 265 / 281 s/traj). Those are an UPPER BOUND: the redo's weak $g^2=1$-$3$ conditions
  the Wilson operator better -> fewer inner-CG iterations -> faster. (Old L2 Nf4/6 < Nf2 because old Nf2 ran
  the serial binary while Nf4/6 used the block force-parallel binary; redo Nf2 = NSTACK1 block = serial-equiv.)

## Per-trajectory (packed 2/GPU) and stream-hours to target
Targets: massless L1 -> 2000 traj/stream, L2 -> 1200; massive L1 -> 120, L2 -> 80.

| group | streams x target | s/traj (est) | stream-hours |
|---|---|---|---|
| massless L1 Nf2 / Nf4 / Nf6 | 3 x 2000 each | 35 / 55 / 70 | ~58 / 92 / 117 (=267) |
| **massless L2 Nf2 / Nf4 / Nf6** | 3 x 1200 each | ~200-400 | **~600-1200** |
| massive L1 Nf2 | 4 x 120 | 35 | ~5 |
| massive L2 Nf2 | 4 x 80 | ~280 | ~25 |

"stream-hours" = sum over streams of (traj target) $\times$ (s/traj). 26 streams total.

## Totals
- **~450-750 A100 GPU-hours** (2 streams/GPU -> GPU-hours $\approx$ stream-hours $/\,2$; add ~10-15% for
  imperfect pair balancing). **~75% of the cost is massless L2** (worst: Nf6).
- **Wall-clock $\approx$ 3-6 days IF all 13 GPUs run uninterrupted** -- set by the heaviest slot (slot 8,
  L2 Nf6 g2.0+g3.0): $1200 \times (200$-$420\,\text{s}) = 3$-$6$ days on one GPU.
- **Massive + all L1 are cheap (~1 day)** and finish first.

## `opp` (preemptible) caveat
`opp` jobs can be preempted mid-run. Trajectories checkpoint (every accepted traj) and jobs `afterany`-chain,
so correctness is unaffected, but the ~6-day CONTINUOUS need on the L2 streams realistically stretches to
**1-2+ weeks** wall-clock under preemption. Cheap slots (massive, L1) barely notice.

## Biggest uncertainty / how to tighten
L2 (esp. Nf6) per-traj is the dominant unknown and is only extrapolated. To nail it: one 30-min `qos=test`
job packing an L2 Nf6 pair -- `QOS=test WALL=00:30:00 WSEC=1800 SLOTS="8" bash run_wrapper_redo_claude.sh 1`
-- read `# HMC : <sec>` from its logs and re-scale the L2 rows.

## Levers if the total is too high
- Trim L2 target 1200 -> 600 (halves the dominant cost).
- Drop / defer Nf6 (the single most expensive Nf).
- Stage: run massive + L1 first (~1 day), commit L2 after seeing real throughput.
