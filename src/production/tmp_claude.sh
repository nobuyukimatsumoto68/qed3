#!/usr/bin/env bash
# tmp_claude.sh -- launch the Nf6 extra-hit run (run_disc_Nf6_hits_claude.sh), detached.
# L1 Nf6 g1.0/g1.5 -> 4 hits, L2 Nf6 g2.0/g3.0 -> 6 hits.  Log: run_disc_Nf6_hits_claude.log.
# GPUS="0 0" bash tmp_claude.sh   would keep GPU1 free.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
echo "### launching run_disc_Nf6_hits_claude.sh  [$(date +%F_%H:%M:%S)] ###" | tee run_disc_Nf6_hits_claude.log
nohup bash run_disc_Nf6_hits_claude.sh >> run_disc_Nf6_hits_claude.log 2>&1 &
echo "### started PID $!  (tail -f run_disc_Nf6_hits_claude.log) ###" | tee -a run_disc_Nf6_hits_claude.log
