# Impl plan -- thorough CONNECTED tower analysis, massless gsq=8, all L=1,2,4 (Nf2/4/6)

## Goal

A dedicated, thorough analysis of the CONNECTED local-current + scalar $Y_{\ell m}$ tower for the **9 massless
gsq=8 ensembles** (Nf 2/4/6 $\times$ L1/L2/L4), at **stride $=2$** (every other config; user 2026-07-14), **nhit $=1$**. Full tower =
vector ($s_1,s_2,s_3$), axial, and scalar ($\sigma_{PS}$, $\sigma_{FS}$), $\ell=0..3$, per $m$. Disc vanishes for the
scalars (flavor/chiral selection at $m=0$), so conn is the full physical scalar correlator.

Definitions are the CURRENT furnished ones (the GW-symmetrized $\sigma_{PS}$ was NOT adopted):
$\sigma_{PS}$ conn $=V_{++}$, $\sigma_{FS}$ conn $=V_{++}+V_{--}^{FS}$ ; vector/axial per `jj_local_ylm_impl_plan_claude.md`.

## Ensembles (survey 2026-07-14; base spacing $=1$ everywhere)

| Nf | L1 avail | L2 avail | L4 avail |
|----|---------:|---------:|---------:|
| 2  | 3426 | 1601 | 168 |
| 4  | 3557 |  920 | 206 |
| 6  | 3163 |  713 |  83 |

At stride $=2$ ($\approx$ 6,900 configs total: L1 $\approx$1713/1779/1582, L2 $\approx$801/460/357, L4 $\approx$84/103/42).
Existing conn (full tower + scalar) at $\sim$140 configs each for L1/L2 is "complete"-gated and will be SKIPPED
where it lands on the stride-2 grid; L4 has none. NOTE: base $=1$ => still autocorrelated at stride 2; analysis
MUST binned-jackknife (bin $\gtrsim \tau_\text{int}$).

## Driver (existing, unchanged)

`jj_local_ylm_scalar_conn_stoch_claude.cu` in DEFAULT (fresh, NON-`--IsScalarOnly`) mode = vector+axial+scalar,
writes `h0/ylm/s{a}`, `h0/ylm_axial/s{a}`, `h0/scalar`, `h0/scalar_fs` + `complete`. Args per ensemble:
`--ens-dir <ens>/ --kmin <first> --stride 2 --kmax <last+1> --nhits 1 --t0 0 --spin-dilution`. Per-hit atomic +
`complete`-gated => resumable. Per-L binaries `-DN_REFINE_CLI={1,2,4}`.

## Files

- NEW `run_conn_gsq8_massless_claude.sh` (HANDOFF, user runs): build L1/L2/L4 binaries (need_build), 9 ensembles,
  every config, cost-balanced MPS pool, auto-detect range, resumable. Logs `conn_massless_L{L}_Nf{nf}_claude.log`.
- NEW `conn_gsq8_massless_analysis_claude.ipynb` (thorough): per-ensemble binned-jackknife; $\sigma_{PS}/\sigma_{FS}$ and
  vector/axial $\ell=0..3$ towers; effective masses; autocorr/bin-size scan; cross-L/cross-Nf comparison. Reads
  `corr_ylm_conn_t00_nhits1_s1/`. Scaffolds on the existing $\sim$140 and grows as the campaign fills in.

## Chunks

1. This plan + `run_conn_gsq8_massless_claude.sh` (measurement handoff). [do first]
2. Analysis notebook scaffold: loaders + binned-jackknife + $\sigma$ towers + effmass on existing data.
3. Thoroughness pass: autocorr/bin scan, vector/axial towers, fits, all-9 cross comparison, as data accumulates.

## Compute / cadence

Full-tower conn per config $\sim$ tens of s (L1) to a few hundred s (L4). $\sim$13.8k configs => multi-day,
multi-GPU; the run script is resumable (kill/restart safe) and MPS-packable. USER runs it; I read logs and build
the analysis. Cost note is surfaced, not hidden.

## Open questions -- RESOLVED (2026-07-14)

- Content = FULL conn tower (scalar + vector + axial).  Density = stride 2 (~6.9k configs).  L4 = INCLUDED.
- nhit=1, t0=0, spin-dilution ON, furnished (non-GW) definitions.  Disc NOT needed (vanishes).