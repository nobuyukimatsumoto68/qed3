#!/bin/bash -l
# tmp_monitor_L4_scc_claude.sh  (_claude handoff, 2026-07-16)
# Throughput monitor for the SCC L4 production run. READ-ONLY: qstat + parse the tee'd run logs and count
# checkpoints. No submit, no rm. Run repeatedly on a login node:  bash tmp_monitor_L4_scc_claude.sh
set -u

cd /projectnb/qfe/nmatsum/qed3/src/production || exit 1

echo "############ SCC L4 monitor  $(date) ############"

echo "==== running / queued L4 jobs (SGE) ===="
qstat -u "$USER" 2>/dev/null | awk 'NR<=2 || /L4[78]0/'
echo

echo "==== per-ensemble throughput (from run_L4_scc_*.log) ===="
printf "%-46s %6s %10s %10s %10s\n" "log" "ntraj" "last_sec" "mean_sec" "last_dH"
for f in $(ls -t run_L4_scc_*.log 2>/dev/null)
do
  ntraj=$(grep -c "# HMC :" "$f" 2>/dev/null)
  last_sec=$(grep "# HMC :" "$f" 2>/dev/null | tail -1 | awk '{print $4}')
  mean_sec=$(grep "# HMC :" "$f" 2>/dev/null | awk '{s+=$4; n++} END{if(n>0) printf "%.1f", s/n; else print "-"}')
  last_dH=$(grep "# dH :" "$f" 2>/dev/null | tail -1 | awk '{print $4}')
  printf "%-46s %6s %10s %10s %10s\n" "$(basename "$f")" "${ntraj:-0}" "${last_sec:--}" "${mean_sec:--}" "${last_dH:--}"
done
echo

echo "==== checkpoint progress vs KMAX (ckpoint_lat.* per output dir) ===="
# output dirs: Nf<Nf>_gsq<g>at0.200000nu01.000000mRe<m>mIm0.000000nt128L4_hb<aux>/
for d in Nf*L4_hb*/
do
  [ -d "$d" ] || continue
  kmax=$(ls "$d"ckpoint_lat.* 2>/dev/null | sed 's/.*\.//' | sort -n | tail -1)
  ncfg=$(ls "$d"ckpoint_lat.* 2>/dev/null | wc -l)
  printf "  %-70s  last_k=%-6s  nckpt=%s\n" "$d" "${kmax:-0}" "$ncfg"
done
echo

echo "==== acceptance summary (running rate, last line per log) ===="
for f in $(ls -t run_L4_scc_*.log 2>/dev/null)
do
  rate=$(grep "rate :" "$f" 2>/dev/null | tail -1 | awk '{print $NF}')
  [ -n "$rate" ] && printf "  %-46s rate=%s\n" "$(basename "$f")" "$rate"
done

echo "############ done ############"
