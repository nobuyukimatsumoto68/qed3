# Ensembles and Analysis Status — gsq=8

All ensembles: gsq=8.0, at=0.2, Nt=128, L=1 (N_REFINE=1, 12 spatial sites).
Gauge configs stored as `ckpoint_lat.K` / `ckpoint_rng.K`.
`ckpoint_lat` files are kept for every trajectory; `ckpoint_rng` files are kept only for the latest trajectory and every 1000th (k_ckpoint_rng=1000), to save disk space while preserving periodic restore points.
Meson observables: HDF5, path `{ell}/{em}/{ab}/{h}/{t0}/real` and `.../imag`, length Nt.
Glueball observables: plain-text `F.K` (16 time-averaged operators) and `F_corr.K` (16x16 correlator matrix).

---

## Quenched (Nf=0)

**Gauge ensemble:**
- Directory: `gsq8.000000at0.200000nt128L1/`
- Configs: 10000 (`ckpoint_lat.10`, `ckpoint_lat.20`, ..., `ckpoint_lat.10000`; step 10)
- Status: complete (long run)

### Meson measurements (quenched, valence mass scan)

Gauge configs read from `gsq8.000000at0.200000nt128L1/`.
Output in `data_gsq8.000000at0.200000nu1{nu1}nt128L1/`.
All runs: nhits=1, dt=128, ellmax=1 (from `tune_nu1_claude.sh` / `tune_nu1_claude_v2.sh`).

| nu1  | # HDF5 files | Traj range       | Notes                  |
|------|-------------|------------------|------------------------|
| 1.25 |  4          | 10, 20, 30, 40   | test/early run only    |
| 1.50 | 100         | 10--1000, step 10 |                       |
| 1.75 | 100         | 10--1000, step 10 |                       |
| 2.00 | 100         | 10--1000, step 10 |                       |
| 2.25 | 100         | 10--1000, step 10 |                       |
| 2.50 | 100         | 10--1000, step 10 |                       |
| 2.75 | 100         | 10--1000, step 10 |                       |

Analysis results:

```
nu1 = 1.50:
nu1 = 1.75:
nu1 = 2.00:
nu1 = 2.25:
nu1 = 2.50:
nu1 = 2.75:
```

Notes:

---

## ==Nf=2, nu0=1.0==

**Gauge ensemble:**
- Directory: `Nf2_gsq8.000000at0.200000nu01.000000nt128L1/`
- Configs: 999 (`ckpoint_lat.1` to `ckpoint_lat.999`)
- Status: kmax=1000 run, stopped at 999

### Glueball measurement

- Directory: `data_Nf2_gsq8.000000at0.200000nu01.000000nt128L1/`
- Files: `F.K` and `F_corr.K`, K=10 to 990, step 10 (99 trajectories)
- Program: `glue2_claude.cu` (gradient flow, Y_{\ell m} projections \ell=0..3, 16 operators)

Analysis results:

```
glueball mass (lightest, \ell=0):
GEVP plateau fit range:
```

Notes: Almost sqrt(2)

### Meson measurement (nu0=nu1=1.0, unitary)

- Directory: `data_Nf2_gsq8.000000at0.200000nu01.000000nu11.000000nt128L1/`
- Files: `meson_corr_v2.K.h5`, K=10 to 990, step 10 (99 HDF5 files)
- Program: `meson_pq_wall_v2_claude.cu`

Analysis results:

```
ab=3 (pseudoscalar) effective mass:
ab=0 (scalar) effective mass:
```

Notes: Almost sqr(2); ***supports no tuning hypothesis***

---

## Nf=2, nu0=1.5

**Gauge ensemble:**
- Directory: `Nf2_gsq8.000000at0.200000nu01.500000nt128L1/`
- Configs: 999 (`ckpoint_lat.1` to `ckpoint_lat.999`)
- Status: kmax=1000 run, stopped at 999

### Meson measurement (nu0=nu1=1.5, unitary)

- Directory: `data_Nf2_gsq8.000000at0.200000nu01.500000nu11.500000nt128L1/`
- Files: `meson_corr_v2.K.h5`, K=10 to 990, step 10 (99 HDF5 files)
- Program: `meson_pq_wall_v2_claude.cu`

Analysis results:

```
ab=3 (pseudoscalar) effective mass:
ab=0 (scalar) effective mass:
```

Notes:

---

## Nf=4, nu0=1.0

**Gauge ensemble:**
- Directory: `Nf4_gsq8.000000at0.200000nu01.000000nt128L1/`
- Configs: 0 (directory exists but is empty)
- Status: not started

Notes:

---

## TODO

- [x] ckpoint_rng pruning implemented in hmc_claude.cu, hmc_fermilab_claude.cu, hmc_fermilab_L2_claude.cu (k_ckpoint_rng=1000)
- [x] Start Nf=4, 6 HMC run
- [ ] Analyze meson correlators: fit effective mass for each nu1 in the quenched scan
- [ ] Compare Nf=2 nu0=1.0 vs nu0=1.5 pseudoscalar mass

Notes:
