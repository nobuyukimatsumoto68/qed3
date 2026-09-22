#!/usr/bin/env python3
# build_free_L2_cache_claude.py
#   Build the free-L2 9-op flavor x geometry sigma^2 cache into an L2-SPECIFIC filename, so it does not
#   collide with / clobber the L1 free cache (the production builder names free caches without encoding L).
#   Reuses the validated matrix_one_config from sigma2_flavorgeom_full_v2.  Free L2 uses the COMPLETE
#   distill_Nv84 basis (single exact config).
#   Run: python3 build_free_L2_cache_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "2")
os.environ.setdefault("NVDIR", "distill_Nv84")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import sigma2_flavorgeom_full_v2_claude as B

OUT = "sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_free_L2_1cfg_d1_claude.npy"


def main():
    print("# building free L2 9-op cache: ENS=%s L=%d NVDIR=%s" % (dc.ENS, dc.L, dc.NVDIR))
    dual = dc.dual_areas_from_mesh()
    C = B.matrix_one_config(dc.KS[0], dual)          # (9,9,DTMAX)
    allC = C[None]                                   # (1,9,9,DTMAX)
    os.makedirs("sigma2_flavor_cache_claude", exist_ok=True)
    np.save(OUT, allC)
    print("# saved -> %s  shape %s" % (OUT, allC.shape))


if __name__ == "__main__":
    main()
