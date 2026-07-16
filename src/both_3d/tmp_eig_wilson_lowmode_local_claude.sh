#!/usr/bin/env bash
# _claude: build (L=1,2,4) + run the Wilson low-mode (lambda_min/lambda_max) distribution scan on the
# massless gsq8 ensembles for Nf in NFS, PACKED on one GPU via CUDA MPS (the scan is latency-bound on
# these small lattices, so the L jobs of one Nf overlap and share SMs). Nf is a RUNTIME arg -> one build
# per L, reused for every Nf. Each (L,Nf) writes its own .dat + .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=eig_wilson_lowmode_claude.log
GPU=0
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4
MPS=1                 # 1 = pack the L jobs of each Nf concurrently via MPS; 0 = serial

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=eig_wilson_lowmode_claude.cu

GSQS=(8 12)           # gsq list to scan (gsq12 = the strong, freeze-prone coupling; deeper low-mode tail).
NFS=(6)               # Nf list. Nf2 gsq8 is already scanning in the current job -> default = just Nf6.
                      #   Use NFS=(2 6) / GSQS=(8 12) for a from-scratch run of all.

# (L stride) -- L4 has few configs, so scan it densely
LEVELS=("1 5" "2 5" "4 1")

{
  echo "######## Wilson low-mode scan  $(date)  GPU=${GPU}  gsq={${GSQS[*]}}  Nf={${NFS[*]}}  MPS=${MPS} ########"

  # ---- build all L once (Nf-independent; nvcc is CPU) ----
  for pair in "${LEVELS[@]}"; do
    set -- $pair
    L="$1"
    echo ""
    echo "==================== BUILD L=${L} ===================="
    if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${L} $SRC $LDFLAGS -o "eig_wilson_lowmode_L${L}.x" ; then
      echo "!!! BUILD FAILED L=${L}"
    fi
  done

  # ---- start the MPS control daemon (must be up BEFORE launch; cannot join retroactively) ----
  if [ "$MPS" -eq 1 ]; then
    if pgrep -f 'nvidia-cuda-mps-control' >/dev/null 2>&1; then
      echo "# MPS control daemon already running"
    else
      echo "# starting MPS control daemon (nvidia-cuda-mps-control -d)"
      nvidia-cuda-mps-control -d || echo "!!! MPS daemon start failed -- jobs will time-slice, not share"
    fi
  fi

  # ---- run: outer loops over gsq and Nf; inner L jobs packed (MPS) or serial ----
  for GSQ in "${GSQS[@]}"; do
  for NF in "${NFS[@]}"; do
    echo ""
    echo "==================== gsq=${GSQ} Nf=${NF} ===================="
    if [ "$MPS" -eq 1 ]; then
      pids=()
      for pair in "${LEVELS[@]}"; do
        set -- $pair
        L="$1"
        strd="$2"
        [ -x "eig_wilson_lowmode_L${L}.x" ] || continue
        jlog="eig_wilson_lowmode_gsq${GSQ}_Nf${NF}_L${L}_claude.log"
        echo "#   gsq=${GSQ} Nf=${NF} L=${L} (stride=${strd}) -> ${jlog}"
        ./"eig_wilson_lowmode_L${L}.x" "$GSQ" "$NF" "$strd" 0 > "$jlog" 2>&1 &
        pids+=($!)
      done
      echo "# waiting for pids: ${pids[*]}"
      wait "${pids[@]}"
      for pair in "${LEVELS[@]}"; do
        set -- $pair
        L="$1"
        echo "==== gsq=${GSQ} Nf=${NF} L=${L} ===="
        cat "eig_wilson_lowmode_gsq${GSQ}_Nf${NF}_L${L}_claude.log" 2>/dev/null
      done
    else
      for pair in "${LEVELS[@]}"; do
        set -- $pair
        L="$1"
        strd="$2"
        [ -x "eig_wilson_lowmode_L${L}.x" ] || continue
        echo "-------------------- RUN gsq=${GSQ} L=${L} Nf=${NF} (stride=${strd}, serial) --------------------"
        ./"eig_wilson_lowmode_L${L}.x" "$GSQ" "$NF" "$strd" 0
      done
    fi
  done
  done

  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
