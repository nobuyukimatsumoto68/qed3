#!/usr/bin/env python3
# t00_ell2_singlet_claude.py
# Chunk C: singlet = C_conn + C_disc from the loops cached by t00_ell2_claude.py (plan: t00_ell2_impl_plan_claude.md).
# Nf=2: O_H = eta^H W xi + h.c. already sums both flavors (W half / W^dag half), so with the driver's definitions
#   C_singlet(dt) = C_conn(dt) + C_disc(dt),  C_disc(dt) = < L_m(t+dt) L_m(t) >_{t,m} - vacuum   (ensemble connected)
# Run: CACHE=t00_ham_cache_claude/ell2_..._claude.npz ELL=2 python3 t00_ell2_singlet_claude.py
import os
import numpy as np

CACHE = os.environ["CACHE"]
ELL = int(os.environ.get("ELL", "2"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
d = np.load(CACHE)
Cc = d["C%d" % ELL]
L = d["L%d" % ELL]
ncfg = L.shape[0]
twin = L.shape[1]
dtmax = Cc.shape[1]
Lbar = L.mean(0)
dL = L - Lbar[None]
Cd = np.zeros((ncfg, dtmax))
for dt in range(dtmax):
    ns = twin - dt
    Cd[:, dt] = (dL[:, dt:dt + ns, :] * dL[:, 0:ns, :]).mean(axis=(1, 2))
nb = ncfg // BINSIZE
bc = np.array([Cc[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
bd = np.array([Cd[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])


def jk(f):
    c = f(bc.mean(0), bd.mean(0))
    s = np.array([f(np.delete(bc, i, 0).mean(0), np.delete(bd, i, 0).mean(0)) for i in range(nb)])
    return c, np.sqrt((nb - 1) * np.nanmean((s - np.nanmean(s, 0)) ** 2, axis=0))


def em(x):
    with np.errstate(all="ignore"):
        return np.log(x[:-1] / x[1:])


d_c, d_e = jk(lambda c, dd: dd)
r_c, r_e = jk(lambda c, dd: dd / c)
ma_c, ma_e = jk(lambda c, dd: em(c))
ms_c, ms_e = jk(lambda c, dd: em(c + dd))
print("# %s  ELL=%d ncfg=%d nb=%d  <L>=%.3e (rms over m,t of cfg-mean)" % (CACHE, ELL, ncfg, nb, np.sqrt((Lbar ** 2).mean())))
print("#  dt |  C_disc(err)              C_disc/C_conn(err)     | m_adj(conn)      m_singlet(conn+disc)")
for dt in range(1, dtmax - 1):
    print("#  %2d | %11.3e(%9.2e)   %9.4f(%8.4f)   | %7.4f(%.4f)  %7.4f(%.4f)"
          % (dt, d_c[dt], d_e[dt], r_c[dt], r_e[dt], ma_c[dt], ma_e[dt], ms_c[dt], ms_e[dt]))
