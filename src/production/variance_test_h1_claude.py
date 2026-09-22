#!/usr/bin/env python3
# variance_test_h1_claude.py  (2026-09-02)
# Noise decomposition for the 2nd stochastic hit (h1) on the connected axial Y_lm l=3 per-m correlator.
# Answers: does adding h1 reduce the variance of the l=3 per-m spectrum, and by how much?
#
# Model (per (m,t)):
#   single-hit estimator per config X0 has variance  sigma_tot^2 = sigma_g^2 (gauge) + sigma_s^2 (stochastic).
#   sigma_s^2 estimated from the two INDEPENDENT hits within a config: E[(X0-X1)^2] = 2 sigma_s^2.
#   sigma_g^2 = sigma_tot^2 - sigma_s^2.
#   n-hit averaged estimator variance = sigma_g^2 + sigma_s^2 / n.
#   predicted 2-hit / 1-hit ERROR ratio r = sqrt( (sigma_g^2 + sigma_s^2/2) / (sigma_g^2 + sigma_s^2) )
#                                          = sqrt( 1 - (1/2) f ),   f = sigma_s^2 / sigma_tot^2 (stochastic fraction).
#   r -> 1/sqrt2 = 0.71 if fully stochastic (f=1); r -> 1 if gauge-dominated (f=0).
# READ-ONLY analysis of existing h5. No writes.

import sys
import glob
import re
import numpy as np
import h5py

ens = sys.argv[1] if len(sys.argv) > 1 else "Nf4_gsq2"
# signal t-window for the summary (source at t0=0; use small-t where the l=3 signal is strong)
TWIN = range(1, 16)

base = glob.glob(f"data_{ens}*L4_hb0.400000-1.000000_v*/corr_ylm_conn_t00_nhits1_s1")
if not base:
    print(f"no conn dir for {ens}")
    sys.exit(1)
cdir = base[0]

MS = [-3, -2, -1, 0, 1, 2, 3]
SPINS = ["s1", "s2"]
VS = ["Vpp", "Vmm"]

def obs(fname):
    # A[m_index, t] = sum over {s1,s2} x {Vpp,Vmm} of the REAL part of the axial l=3 (l,m) correlator.
    # (a fixed linear combination; the stochastic-fraction result is robust to the exact combo.)
    out = np.zeros((len(MS), 128))
    with h5py.File(fname, "r") as f:
        g = f["h0/ylm_axial"]
        for mi, m in enumerate(MS):
            acc = np.zeros(128)
            for s in SPINS:
                for V in VS:
                    acc += np.array(g[f"{s}/l3/m{m}/{V}/real"])
            out[mi] = acc
    return out

def kof(p):
    return int(re.search(r"corr\.(\d+)\.h", p).group(1))

h0_files = {kof(p): p for p in glob.glob(f"{cdir}/corr.*.h0.h5")}
h1_files = {kof(p): p for p in glob.glob(f"{cdir}/corr.*.h1.h5")}
both = sorted(set(h0_files) & set(h1_files))
allh0 = sorted(h0_files)
print(f"# ensemble {ens}: full h0 configs = {len(allh0)}, both-hit (h0&h1) configs = {len(both)}")
if len(both) < 5:
    print("# too few both-hit configs for a stable sigma_s^2 estimate")
    sys.exit(0)

# sigma_s^2 (per (m,t)) from the both-hit subset: mean_c (X0-X1)^2 / 2
d2 = np.zeros((len(MS), 128))
for k in both:
    a0 = obs(h0_files[k])
    a1 = obs(h1_files[k])
    d2 += (a0 - a1) ** 2
sig_s2 = 0.5 * d2 / len(both)

# sigma_tot^2 (single-hit config variance) from the FULL h0 set
stack = np.stack([obs(h0_files[k]) for k in allh0], axis=0)  # (Ncfg, m, t)
sig_tot2 = stack.var(axis=0, ddof=1)

sig_g2 = np.clip(sig_tot2 - sig_s2, 0.0, None)
with np.errstate(divide="ignore", invalid="ignore"):
    f_stoch = np.where(sig_tot2 > 0, sig_s2 / sig_tot2, np.nan)
    ratio = np.sqrt(np.clip((sig_g2 + 0.5 * sig_s2) / (sig_g2 + sig_s2), 0, 1))

tw = list(TWIN)
print(f"# t-window for summary = {tw[0]}..{tw[-1]}")
print(f"# {'m':>3} {'stoch_frac_f':>13} {'2hit/1hit_err':>14} {'eff_gain%':>10}")
for mi, m in enumerate(MS):
    ff = np.nanmean(f_stoch[mi, tw])
    rr = np.nanmean(ratio[mi, tw])
    gain = 100.0 * (1.0 - rr)
    print(f"  {m:>3} {ff:>13.3f} {rr:>14.3f} {gain:>9.1f}%")
allf = np.nanmean(f_stoch[:, tw])
allr = np.nanmean(ratio[:, tw])
print(f"# l=3 ALL-m: stoch_frac={allf:.3f}  2hit/1hit_err={allr:.3f}  (err reduction {100*(1-allr):.1f}%)")
print("# interp: f~1 & ratio~0.71 => stochastic-limited, h1 big win;  f~0 & ratio~1.0 => gauge-limited, h1 useless.")
