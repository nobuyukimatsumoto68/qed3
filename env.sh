module load cuda/12.8
# module load python3
# module load llvm
module load gcc/13.2.0

# --- unified qed3 build environment (single source of truth; source this file, nothing else) ---
# Repo root = the directory containing this env.sh (robust to relocation).
export QED3_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
# The HMC / production CUDA build needs ONLY Eigen (installed in the repo by set.sh: qfe_mod/include/Eigen)
# on top of CUDA. HDF5 / HighFive / GSL are for the measurement programs (glue, jj) and are deliberately
# NOT loaded here -- out of scope for the HMC ensemble generation.
export QED3_INC="-I${QED3_ROOT}/qfe_mod/include"
