#!/usr/bin/env python3
# distill_contract_claude.py
#
# CHUNK 3 of exact distillation (distillation_impl_plan_claude.md): contract the perambulator into the
# scalar one-point loops and VALIDATE against the stochastic loops on the SAME configs (gold test T3a).
#
# Vertex (matched EXACTLY to valence_claude.h mult_Ylm_real(0,0) + accumulate_loop_raw):
#   W_0 = diag_x( dual_areas[x] * Y00 ) (x) 1_spin ,  Y00 = 1/sqrt(4 pi) ,  NO 1/N_SITES.
#   Phi(t) = V(t)^dag W_0 V(t)   (Nv x Nv).
# Distilled loops from FORWARD perambulators (spin folded into the mode index):
#   DS(t)      = Tr[Phi tau(t,t)]                      (= stochastic DS = J_1)
#   DS_1mD(t)  = Tr[Phi tau(t,t)] - Tr[Phi]            (= tr[W_0 (1-D_ov) D_ov^{-1}], since (1-D_ov)G = G-1)
#   Dp_est(t)  = Tr[Phi tau(t,t)^2]                    (matches the STOCHASTIC estimator: bare middle Pi_t)
#   Dp_phys(t) = Tr[Phi tau(t,t) Phi tau(t,t)]         (PHYSICAL connected sigma\sigma: TWO area vertices)
#   Dp1mD(t)   = Tr[Phi tau_gw(t,t)^2]                 (FS extended, estimator form)
#
# NOTE (vertex-bookkeeping pin): the stochastic D'_S estimator (jj_sigma_loops) uses a BARE timeslice
# projector in the middle -> it computes Dp_est = Tr[Phi tau^2], NOT the physical Dp_phys = Tr[Phi tau Phi
# tau].  At L1 the icosahedron is vertex-transitive so dual_areas are UNIFORM and Phi = w00 * 1, giving
# Dp_phys = w00 * Dp_est exactly -- so they differ only by one vertex weight w00.  This script prints both
# and the w00 factor so the physical operator can be pinned vs v3-4 Ch.5.
#
# T3a: distilled DS/DS_1mD (EXACT, per config) vs the stochastic corr_sigma_loops (noisy, hit-averaged) on
# the same configs -- must agree within the stochastic scatter and tighten it.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import glob
import re
import numpy as np
import h5py
import math

GEOM = "../../geometry/data/"
L = int(os.environ.get("LREF", "1"))          # refinement level; LREF=2 for the L2 free test
NS = 2
Y00 = 1.0 / math.sqrt(4.0 * math.pi)

ENS = os.environ.get("ENS", "Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
NVDIR = os.environ.get("NVDIR", "distill_Nv24")   # NVDIR=distill_Nv84 for the L2 free (complete) basis
PERAM_DIR = "data_" + ENS + "/" + NVDIR + "/"
LOOP_DIR = "data_" + ENS + "_vmRe0.000000vmIm0.000000/corr_sigma_loops_tb2_nhits2/"


def discover_ks():
    # all configs with a perambulator whose /peram/tau is present, sorted
    ks = []
    for f in glob.glob(PERAM_DIR + "peram.*.h5"):
        k = int(re.search(r"peram\.(\d+)\.h5", f).group(1))
        try:
            with h5py.File(f, "r") as h:
                if "peram" in h and "tau" in h["peram"]:
                    ks.append(k)
        except Exception:
            pass
    return sorted(ks)


KS = discover_ks()


def load_vec3(path):
    out = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if len(s) < 3:
                continue
            out.append([float(s[0]), float(s[1]), float(s[2])])
    return np.array(out)


def load_intlists(path):
    out = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            out.append([int(v) for v in s])
    return out


def sph_tri_area(p, x, y):
    # spherical triangle area (excess) for unit vectors p,x,y -- matches set_dual_areas (IS_FLAT off).
    a = math.acos(max(-1.0, min(1.0, float(np.dot(x, p)))))
    b = math.acos(max(-1.0, min(1.0, float(np.dot(y, p)))))
    c = math.acos(max(-1.0, min(1.0, float(np.dot(x, y)))))
    s = 0.5 * (a + b + c)
    t = math.tan(0.5 * s) * math.tan(0.5 * (s - a)) * math.tan(0.5 * (s - b)) * math.tan(0.5 * (s - c))
    return 4.0 * math.atan(math.sqrt(max(t, 0.0)))


def dual_areas_from_mesh():
    sites = load_vec3(GEOM + "pts_n%d.dat" % L)
    dual_sites = load_vec3(GEOM + "pts_dual_n%d.dat" % L)
    dual_faces = load_intlists(GEOM + "face_dual_n%d.dat" % L)
    n = sites.shape[0]
    area = np.zeros(n)
    for ip in range(n):
        p = sites[ip]
        face = dual_faces[ip]
        acc = 0.0
        m = len(face)
        for i in range(m):
            ix = face[i]
            iy = face[(i + 1) % m]
            acc += sph_tri_area(p, dual_sites[ix], dual_sites[iy])
        area[ip] = acc
    return area


def _read_peram_raw(f):
    # Reads either layout and NORMALIZES tau/tau_gw to a LEADING source-window axis:
    #   nsrc==1 (old _v1 / --nsrc 1): peram/tau = (twin,twin,Nv,Nv) -> promoted to (1,twin,twin,Nv,Nv)
    #   nsrc>1  (--nsrc 2, L2):       peram/tau = (nsrc,twin,twin,Nv,Nv) verbatim
    # Returns (V, tau, taugw, tsrc_list, twin) with tau shape (nsrc, twin, twin, Nv, Nv).
    V = f["V/real"][:] + 1j * f["V/imag"][:]       # (Nt, Nv, 2Ns) -- all-t, SHARED across windows
    tau = f["peram/tau/real"][:] + 1j * f["peram/tau/imag"][:]
    taugw = f["peram/tau_gw/real"][:] + 1j * f["peram/tau_gw/imag"][:]
    twin = int(f["meta/twin"][0])
    if tau.ndim == 4:                              # nsrc==1 -> add the leading window axis
        tau = tau[None]
        taugw = taugw[None]
    if "meta/tsrc_list" in f:                       # per-window base source timeslices
        tsrc_list = [int(x) for x in f["meta/tsrc_list"][:]]
    else:                                           # old _v1 files: single window at meta/tsrc0
        tsrc_list = [int(f["meta/tsrc0"][0])]
    return V, tau, taugw, tsrc_list, twin


def load_peram(k):
    # BACKWARD-COMPAT: returns the FIRST source window (V, tau, taugw, tsrc0, twin), exactly as before.
    # Works unchanged on old _v1 / --nsrc 1 files AND on multi-window --nsrc 2 files (window 0).
    with h5py.File(PERAM_DIR + "peram.%d.h5" % k, "r") as f:
        V, tau, taugw, tsrc_list, twin = _read_peram_raw(f)
    return V, tau[0], taugw[0], tsrc_list[0], twin


def load_peram_windows(k):
    # ALL source windows -> (V, windows) with windows = [(tsrc0_s, tau_s, taugw_s), ...], one per window s.
    # tau_s/taugw_s are (twin,twin,Nv,Nv); V is shared (all-t).  For nsrc==1 files this is a 1-element list,
    # so callers can treat every file uniformly.  Combine windows as INDEPENDENT source sets (~2x configs @ L2).
    with h5py.File(PERAM_DIR + "peram.%d.h5" % k, "r") as f:
        V, tau, taugw, tsrc_list, twin = _read_peram_raw(f)
    return V, [(tsrc_list[s], tau[s], taugw[s]) for s in range(tau.shape[0])]


def load_loop(k, key):
    # average the two hits of the stochastic loop key ("DS","DS_1mD","Dp","Dp_1mD") at l0/m0.
    acc = None
    nh = 0
    for h in (0, 1):
        p = LOOP_DIR + "corr.%d.h%d.h5" % (k, h)
        if not os.path.exists(p):
            continue
        with h5py.File(p, "r") as f:
            base = "h0/sigma_loops/%s/l0/m0/J" % key
            re = np.array(f[base + "/real"])
            im = np.array(f[base + "/imag"])
        v = re + 1j * im
        acc = v if acc is None else acc + v
        nh += 1
    return acc / nh if nh else None


def main():
    dual = dual_areas_from_mesh()
    n_sites = dual.shape[0]
    print("# dual_areas: N_sites=%d  sum=%.10f (4pi=%.10f)  min=%.6e max=%.6e  uniform_dev=%.3e"
          % (n_sites, dual.sum(), 4.0 * math.pi, dual.min(), dual.max(),
             (dual.max() - dual.min()) / dual.mean()))
    # spinor-space vertex weight w00_spinor[j] = dual_areas[j//2] * Y00  (j = NS*x + s)
    w00 = np.repeat(dual, NS) * Y00                 # (2Ns,)
    w00_uniform = dual.mean() * Y00
    print("# w00 (uniform value, for the Dp_phys/Dp_est factor) = %.10f" % w00_uniform)

    for k in KS:
        V, tau, taugw, tsrc0, twin = load_peram(k)
        Nv = V.shape[1]
        # distilled loops over the window
        DS = np.zeros(twin, complex)
        DS1mD = np.zeros(twin, complex)
        Dp_est = np.zeros(twin, complex)
        Dp_phys = np.zeros(twin, complex)
        Dp1mD = np.zeros(twin, complex)
        for a in range(twin):
            t = tsrc0 + a
            Vt = V[t].T                             # (2Ns, Nv), column k = mode k
            Phi = (Vt.conj().T) @ (w00[:, None] * Vt)   # (Nv, Nv)
            tt = tau[a, a]                          # tau(t,t)  (Nv,Nv)
            gg = taugw[a, a]                        # tau_gw(t,t)
            trPhi = np.trace(Phi)
            DS[a] = np.trace(Phi @ tt)
            DS1mD[a] = DS[a] - trPhi
            Dp_est[a] = np.trace(Phi @ tt @ tt)
            Dp_phys[a] = np.trace(Phi @ tt @ Phi @ tt)
            Dp1mD[a] = np.trace(Phi @ gg @ gg)
        # stochastic loops (hit-averaged), restricted to the window
        sDS = load_loop(k, "DS")
        if sDS is None:
            print("\n=== k=%d : no stochastic loop file (T3a skipped) ===" % k)
            continue
        sDS1mD = load_loop(k, "DS_1mD")
        sDp = load_loop(k, "Dp")
        sDp1mD = load_loop(k, "Dp_1mD")
        w = slice(tsrc0, tsrc0 + twin)
        print("\n=== k=%d  (Nv=%d, tsrc0=%d, twin=%d) ===" % (k, Nv, tsrc0, twin))
        # T3a: distilled vs stochastic DS across the window
        rDS = DS / sDS[w]
        rDS1 = DS1mD / sDS1mD[w]
        print("  [T3a DS]     Re DS(t=0..4) distilled = %s" % np.array2string(DS[:5].real, precision=5))
        print("               Re DS(t=0..4) stochastic= %s" % np.array2string(sDS[w][:5].real, precision=5))
        print("               DS  distilled/stochastic (Re): mean=%.4f std=%.4f over window"
              % (rDS.real.mean(), rDS.real.std()))
        print("               DS_1mD distilled/stochastic (Re): mean=%.4f std=%.4f"
              % (rDS1.real.mean(), rDS1.real.std()))
        # Dp estimator match + physical vs estimator factor
        rDp = Dp_est / sDp[w]
        print("  [Dp est]     Dp_est distilled/stochastic (Re): mean=%.4f std=%.4f" % (rDp.real.mean(), rDp.real.std()))
        rphys = Dp_phys / Dp_est
        print("  [Dp pin]     Dp_phys/Dp_est (Re): mean=%.6f  (expect w00=%.6f at L1 uniform)"
              % (rphys.real.mean(), w00_uniform))
        rDp1 = Dp1mD / sDp1mD[w]
        print("  [Dp_1mD est] Dp_1mD distilled/stochastic (Re): mean=%.4f std=%.4f" % (rDp1.real.mean(), rDp1.real.std()))


W10 = np.array([4, 2, 4, 4, 2, 1, 4, 1, 1, 1])
DIAG_LABELS = ["A(-S_S)", "B(-T_S)", "C", "D", "E(C_S^2)", "F", "G", "H", "I", "J"]


def compute_diags(Phi, legs, twin):
    # 10-diagram sum for ONE leg object (legs[a_snk,a_src] = tau, or -taugw for the furnished FS legs).
    # Source fixed at a_src=0 (t=tsrc0), sink a_snk=dt.  Returns diags (10,twin) and CS (twin).
    Phi0 = Phi[0]
    t00 = legs[0, 0]
    DS0 = np.trace(Phi0 @ t00)
    DpS0 = np.trace(Phi0 @ t00 @ Phi0 @ t00)
    diags = np.zeros((10, twin), complex)
    CS = np.zeros(twin, complex)
    for dt in range(twin):
        Phit = Phi[dt]
        ttt = legs[dt, dt]
        t0t = legs[0, dt]                       # (0<-t) leg
        tt0 = legs[dt, 0]                       # (t<-0) leg
        DSt = np.trace(Phit @ ttt)
        DpSt = np.trace(Phit @ ttt @ Phit @ ttt)
        M = Phi0 @ t0t @ Phit @ tt0             # C_S transfer matrix
        CS[dt] = np.trace(M)
        TS = np.trace(M @ M)
        VS_0t = np.trace(Phi0 @ t00 @ M)
        SS_0t = np.trace(Phi0 @ t00 @ Phi0 @ t0t @ Phit @ ttt @ Phit @ tt0)
        VS_t0 = np.trace(Phit @ ttt @ Phit @ tt0 @ Phi0 @ t0t)
        diags[0, dt] = -SS_0t
        diags[1, dt] = -TS
        diags[2, dt] = DSt * VS_0t
        diags[3, dt] = DS0 * VS_t0
        diags[4, dt] = CS[dt] ** 2
        diags[5, dt] = DpS0 * DpSt
        diags[6, dt] = -DS0 * DSt * CS[dt]
        diags[7, dt] = -DS0 ** 2 * DpSt
        diags[8, dt] = -DSt ** 2 * DpS0
        diags[9, dt] = DS0 ** 2 * DSt ** 2
    return diags, CS


def four_point():
    # PS.PS = 2 G10[tau];  FS.FS = G10[tau] + G10[-taugw].  Per config, then average.
    dual = dual_areas_from_mesh()
    w00 = np.repeat(dual, NS) * Y00
    PSPS_all = []
    FSFS_all = []
    CSp_all = []
    CSf_all = []
    ncfg = 0
    for k in KS:
        V, tau, taugw, tsrc0, twin = load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        diags_ps, CSp = compute_diags(Phi, tau, twin)            # S_4 (plain legs) = PS S~_4 = PS S_4
        diags_fs, CSf = compute_diags(Phi, -taugw, twin)         # FS S~_4 (furnished legs)
        S4_ps = np.tensordot(W10, diags_ps, axes=(0, 0))
        S4_fs = np.tensordot(W10, diags_fs, axes=(0, 0))
        PSPS_all.append(2.0 * S4_ps)
        FSFS_all.append(S4_ps + S4_fs)
        CSp_all.append(CSp)
        CSf_all.append(CSf)
        ncfg += 1
    PSPS = np.mean(PSPS_all, axis=0)
    FSFS = np.mean(FSFS_all, axis=0)
    CSp = np.mean(CSp_all, axis=0)
    CSf = np.mean(CSf_all, axis=0)

    print("\n=== two-meson four-point (10 diagrams, %d configs) ===" % ncfg)
    print("  PS.PS = 2 G10[tau] ;  FS.FS = G10[tau] + G10[-tau_gw]")
    for name, X in (("PS.PS", PSPS), ("FS.FS", FSFS)):
        mi = np.max(np.abs(X.imag) / (np.abs(X.real) + 1e-300))
        # connected (subtract the dt->inf plateau estimated from the last window slices)
        plateau = np.mean(X.real[-4:])
        conn = X.real - plateau
        print("  [%s] Re(dt=0..8) = %s" % (name, np.array2string(X[:9].real, precision=4)))
        print("        max|Im/Re|=%.2e   plateau~%.4f   connected(dt=1,2,4)= %.4e %.4e %.4e"
              % (mi, plateau, conn[1], conn[2], conn[4]))
    # single-meson masses from C_S (PS) and C_S^FS, log-effmass over the L1 scalar plateau
    with np.errstate(all="ignore"):
        mp = np.log(CSp[:-1].real / CSp[1:].real)
        mf = np.log(CSf[:-1].real / CSf[1:].real)
    lo, hi = 8, 24
    print("  [C_S PS ] log-effmass dt[%d,%d] = %.4f +- %.4f" % (lo, hi, np.nanmean(mp[lo:hi]), np.nanstd(mp[lo:hi])))
    print("  [C_S FS ] log-effmass dt[%d,%d] = %.4f +- %.4f" % (lo, hi, np.nanmean(mf[lo:hi]), np.nanstd(mf[lo:hi])))
    print("  NOTE: %d configs = machinery/validation only; the coupled GEVP + Delta_- need the L1 sweep." % ncfg)


if __name__ == "__main__":
    main()
    four_point()
