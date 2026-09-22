# Fused multi-source perambulator solve (nsrc>1 at nsrc=1 solve cost) -- impl plan

## Goal / physics
The perambulator solve $\psi_l = D_\text{ov}^{-1} w_l(t)$ produces the FULL-volume propagator. Two source
timeslices that are well separated on the $N_t=128$ torus (e.g. $t_0=0$ and $t_0+N_t/2=64$) do not interfere
in each other's window: by linearity $D_\text{ov}^{-1}(w_l(0)+w_l(64)) = \psi_0 + \psi_{64}$, and in window 0
($t'\in[0,\text{twin})$) the $\psi_{64}$ tail is suppressed by $e^{-m_\text{gap}\,\Delta}$, $\Delta\ge 44$ for
the analysis range $t<20$. So ONE block solve of a COMBINED source yields BOTH windows' perambulators ->
nsrc=2 at the solve cost of nsrc=1 (the L2 wall-time win). Solves stay = twin block-solves (not nsrc*twin).

Idea + folding note from NM: analysis uses $t<20$ so the far-source tail ($\Delta\ge44$) is negligible; the
same solve also carries backward propagation (time-fold) -- deferred to Phase 2.

## Key facts (verified in distill_peram_mrhs_claude.cu)
- Current build (`:532` `for s<nsrc`, `:546` `for a<twin`, `:556` one `blk_Dsq.solve_sq_from_cpu` per source $t$)
  RE-SOLVES per window -> nsrc*twin solves.  `:560` contracts sinks only in $[src0, src0+twin)$.
- `embed_mode_source(src,V_t,l,t)` (`:189`) memsets then writes mode $l$ at ONE timeslice $t$.
- Solve does $(D^\dagger D)^{-1}D^\dagger w = D^{-1}w$ (`:554-556`); linear in the source, so a summed source is fine.
- Output h5 layout is `(nsrc,twin,twin,Nv,Nv)` -- UNCHANGED by fusion (same Tau_all filled, just cheaper).

## STATUS (2026-09-11)
Phase 1 IMPLEMENTED in `distill_peram_mrhs_claude.cu` (backward-compatible: `fused=false` runs the original
per-window path verbatim, guarded by `if(!fused)`; the running L1 nsrc2 job is unaffected). Added:
`embed_mode_source_fused`, `--fused-sources` flag, the fused pre-loop solve, and a `[fused-check]` that reports
the cross-source contamination (full-window and analysis-range `ap,a<20`) instead of asserting. Build+validate
handoff = `tmp_claude.sh` -> `fused_validate_claude.log` (builds a SEPARATE binary
`distill_peram_mrhs_fused_claude.o`, runs one L1 config nsrc=2 with/without `--fused-sources`, compares).
Phase 2 (backward fold) still deferred.

## Phase 1 -- fused-source solve (this change)
Add `--fused-sources` (bool, default off; getopt long only).  When ON and nsrc>1:
- ONE loop over offset $a\in[0,\text{twin})$.  Build block col $l$ = $D_\text{ov}^\dagger\big(\sum_s w_l(src0_s+a)\big)$
  via a new `embed_mode_source_fused` (memset once, add mode $l$ at every window's $src0_s+a$ from its own
  basis $V(src0_s+a)$).  ONE `blk_Dsq.solve_sq_from_cpu` -> $\psi = \sum_s\psi_{src0_s+a}$.
- For each window $s$, extract sinks $tp=src0_s+ap$, $ap\in[0,\text{twin})$: `Tau_all[s][ap][a].col(l)=V(tp)^dag psi|_{tp}`
  (and Taup with $(1-D^\dagger)\psi$).  The other windows' tails are the (tiny) contamination.
- Solves per config = `twin` (vs `nsrc*twin`).  When OFF -> the EXISTING per-window path runs verbatim (preserved).

Files: EDIT `distill_peram_mrhs_claude.cu` -- add `embed_mode_source_fused`; add `--fused-sources` parse + a
`fused` bool; add the fused branch alongside the current per-window loop (original kept, NOT deleted).

## Validation (one config, handoff `tmp_claude.sh` -> `*_claude.log`)
1. Build the binary.
2. Run config k with `--nsrc 2 --twin 32` WITHOUT `--fused-sources` -> peram.k.h5 (out-suffix `_sep`).
3. Run the SAME k WITH `--fused-sources` -> peram.k.h5 (out-suffix `_fused`).
4. Python compare (analysis-ok to run): for each window s, `max|tau_fused - tau_sep|` and the same for tau_gw,
   restricted to the analysis sink range $ap<20$ AND over the full window.  Expect: max-dev over $ap<20$ at the
   contamination floor ($\sim e^{-m_\text{gap}\cdot 44}$, tiny); full-window dev may be larger near the window
   edge closest to the other source -- report both so NM sees the tail.
5. Speedup: wall-time of the fused run vs the separate nsrc=2 run (expect ~2x fewer solves -> ~2x).

## Phase 2 -- backward-window / time-fold extraction (DEFERRED, separate change + validation)
Extract sinks $t' < src0$ (the $-dt$ half) from the same $\psi$ and fold.  Needs its own window-placement +
forward/backward-consistency validation; propose after Phase 1 validates.

## Open questions
- Second source offset: fixed $N_t/2=64$ (from `--nsrc 2` default `tsrc_list=[0,64]`) -- keep using `tsrc_list`
  so `--tsrc-list a,b` still controls placement; fusion just solves them together.  OK?
- Contamination acceptance threshold: below the perambulator solver `tol` (massless run tol=1e-8) over $ap<20$?

Refs: distillation Peardon 0905.2160; mrhs block CG Jegerlehner hep-lat/9612014; multi-source dilution
(well-separated supports) standard; parent `distill_peram_mrhs_impl_plan_claude.md`.
