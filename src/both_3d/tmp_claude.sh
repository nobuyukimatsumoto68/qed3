#!/usr/bin/env bash
# Recompile (current source -- includes the CSR-bucketing sparse_dirac_claude.h) + run the FREE-FIELD
# (U=1) STOCHASTIC exact-conserved-current K correlator at L=1, to validate it against the deterministic
# ground truth corr_deter_exactsum_L1 (same Eq.(4.29)/(4.32)/(4.36) observable -- see
# stoch_jj_switch_handoff_claude.md, Sec.1).
#
#   - free mode  = run_exactjj_free omits --ens-dir (U=1).
#   - NHITS=64   (handoff: conn tp ~6-7 sigma at 64 hit-samples).  Per-hit output + resume:
#                 it SKIPS the 12 already-complete hits in corr_nt02_nhits64/ (those are physics-
#                 identical: only the byte-identical CSR include changed since they were written),
#                 and fills h12..h63.  Partial results are usable at any time.
#   - n_t0=2.
#   - GPU: arg 1 (default 1; L=2 exactsum was stopped so GPU1 is free).  exactsum L=1 still occupies
#                 GPU0, so keep this on GPU1 to avoid compute contention.
#
# Output: data_free_vmRe0.000000vmIm0.000000/corr_nt02_nhits64/corr.0.h{0..63}.h5
# Then:   jj_free_stoch_validate_claude.ipynb  (CFT checks: G_s/G_t -> -1 with NO /(D-1); ylm rates 2,3;
#         ratio 2.4; g1/Gt -> 1/3; g0 -> 0), overlaying corr_deter_exactsum_L1 once it completes.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

GPU="${1:-1}"
NHITS="${2:-64}"
LOG=jj_free_stoch_L1_claude.log

echo "### recompile + run FREE stochastic exact-K jj  L=1  nhits=${NHITS}  GPU=${GPU}  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
bash run_exactjj_free_claude.sh "1" "${NHITS}" "${GPU}" 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}
echo "### run_exactjj_free exit status = ${st}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
exit "$st"
