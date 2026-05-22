# Summary of CUDA Programs and Jupyter Notebooks

This file documents the active source files in `/mnt/barracuda22/qed3/qed3/src/both_3d/`.
Compile-time lattice parameters (N_REFINE, Nt, etc.) are set as `constexpr` inside each file.
Runtime coupling and mass parameters are passed on the command line.

---

## CUDA Programs (_claude.cu)

### disc_claude.cu

Computes free-field meson two-point functions using the overlap Dirac operator on an S2 simplicial lattice (N_REFINE=1, Nt=96).
Two independent Z2 stochastic sources are inverted with the overlap operator, and the resulting propagators are contracted with each of the four gamma matrices (identity and sigma_1,2,3) to produce meson correlators.
Output files `meson_{igam}_corr.free` are written to `data_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nt{Nt}L{L}/`.
The gauge-config loop is commented out; as compiled the program runs only on a trivial (free) gauge background.

CLI: `./a.out [gsq] [Nf] [nu0]`
- `gsq`  Wilson coupling squared (default: 1.0)
- `Nf`   number of fermion flavors (default: 2)
- `nu0`  mass parameter \nu_0 (default: 1.0)
- `-h`   print help and exit

Notes: To be edited; for disconnected diagram calculations

---

### exact_current_claude.cu

Development/scratch program for computing an exact lattice current matrix element using the overlap Dirac operator (N_REFINE=1, Nt=128).
It builds a COO sparse matrix from the Wilson-Dirac operator at a fixed link and applies it to a random source vector, which is the first step toward an exact current insertion calculation.
As written, the program tests the COO matrix-vector product then exits; no output files are produced.
This is an incomplete tool rather than a production measurement.

CLI: `./a.out [options]`
- `--gsq <gsq>`       Wilson coupling squared (default: 2.0)
- `--Nf <Nf>`         number of fermion flavors (default: 2)
- `--nu0 <nu0>`       valence mass \nu_0 (default: 1.0)
- `--nu1 <nu1>`       valence mass \nu_1 (default: 1.0)
- `--nhits <nhits>`   number of stochastic hits (default: 1)
- `--dt <dt>`         time slice separation (default: Nt/2)
- `--ellmax <ellmax>` max angular momentum \ell (default: 2)
- `-h, --help`        print help and exit

Notes: To be edited; for disconnected diagram calculations

---

### gauge_claude.cu

Pure-gauge HMC simulation (no dynamical fermions) on the S2 simplicial lattice (N_REFINE=1, Nt=128, nsteps=20, gsq=8.0, at=0.2 hardcoded).
Accepts an optional directory prefix to select which run to continue; gauge configurations and RNG state are checkpointed every 10 trajectories to `gsq{gsq}at{at}nt{Nt}L{L}/ckpoint_lat.K` and `ckpoint_rng.K`.
The run auto-resumes from the latest checkpoint.

CLI: `./a.out [prefix]`
- `prefix`  directory prefix for reading/writing gauge configs (default: `""`)
- `-h`      print help and exit

Notes:

---

### ==glue2_claude.cu==

Glueball correlator measurement using gradient flow on saved gauge configurations (N_REFINE=1, Nt=128).
For each checkpointed configuration, applies gradient flow (step size 1.0, 100 steps) and projects the flowed plaquette onto real spherical harmonics Y_{\ell m} for \ell=0,1,2,3 (all m), computing 16 operators per time slice.
Two output files are written per trajectory k: `data_{...}/F.K` (time-averaged operator vector, 16 floats) and `F_corr.K` (full 16x16 temporal correlator matrix at all separations dt).
If Nf=0, reads from a pure-gauge directory; otherwise reads from an Nf-flavor directory.

CLI: `./a.out [gsq] [Nf] [nu0]`
- `gsq`  Wilson coupling squared (default: 8.0)
- `Nf`   number of fermion flavors; if 0, reads pure-gauge directory (default: 2)
- `nu0`  mass parameter \nu_0 used in directory name (default: 1.0)
- `-h`   print help and exit

Notes: See run_glue_example.sh

---

### glue_pg_claude.cu

Pure-gauge glueball correlator measurement with gradient flow (N_REFINE=1, Nt=96).
Reads checkpointed pure-gauge configurations from `gsq{gsq}at{at}nt{Nt}L{L}/ckpoint_lat.K`, applies gradient flow, and computes three observables per time slice: the raw plaquette average, the flow-plaquette average, and the chair-angle average.
Per-trajectory outputs: `data_gsq{...}/F.K` (3 time-averaged observables) and `F_corr.K` (3x3 temporal correlator matrix at all separations; the separation index dt is printed as the first column).

CLI: `./a.out [gsq]`
- `gsq`  Wilson coupling squared (default: 12.0)
- `-h`   print help and exit

Notes:

---

### ==hmc_claude.cu==

Full HMC simulation with dynamical overlap fermions on the S2 simplicial lattice (N_REFINE=1, Nt=128, Zolotarev order 11).
Uses the MinimumNorm2 (Omelyan-Mryglod-Folk) integrator with tmax=1.9 and nsteps depending on Nf (4 for Nf=2, 5 for Nf=4 or Nf=6).
Gauge configurations and RNG state are saved every trajectory to `Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nt{Nt}L{L}/ckpoint_lat.K`; the run auto-resumes from the latest checkpoint (kmax=1000).
Per-trajectory stdout reports dH and acceptance rate.

CLI: `./a.out [gsq] [Nf] [nu0]`
- `gsq`  Wilson coupling squared (default: 8.0)
- `Nf`   number of fermion flavors, must be even (default: 2)
- `nu0`  sea mass parameter \nu_0 (default: 1.0)
- `-h`   print help and exit

Notes: **Call in local from run_wrapper_nf_v2.sh** ckpoint_rng files are deleted after each trajectory except every 1000th (k_ckpoint_rng=1000); ckpoint_lat files are all kept.

---

### hmc_fermilab_claude.cu

Full HMC simulation with dynamical overlap fermions, identical in structure to `hmc_claude.cu` but compiled with NPARALLEL_GAUGE=4 for improved parallel throughput and Zolotarev order 21 for higher precision (N_REFINE=1, Nt=128).
kmax is set to 20 (intended as a short benchmark/test run rather than a production run).
All other parameters (tmax=1.9, MinimumNorm2 integrator, nsteps logic by Nf, checkpointing scheme) are identical to `hmc_claude.cu`.

CLI: `./a.out [gsq] [Nf] [nu0]`
- `gsq`  Wilson coupling squared (default: 8.0)
- `Nf`   number of fermion flavors, must be even (default: 2)
- `nu0`  sea mass parameter \nu_0 (default: 1.0)
- `-h`   print help and exit

Notes: ckpoint_rng files are deleted after each trajectory except every 1000th (k_ckpoint_rng=1000); ckpoint_lat files are all kept.

---

### hmc_fermilab_L2_claude.cu

HMC simulation identical to `hmc_fermilab_claude.cu` but for the finer L=2 lattice (N_REFINE=2, 42 spatial sites) with NPARALLEL_DUPDATE=4 for 4-stream parallel overlap operator updates.
Default coupling is gsq=12.0 (stronger coupling appropriate for the smaller lattice spacing at N_REFINE=2).
kmax=20, so this is also a short benchmark run; Zolotarev order is 21, tmax=1.9.

CLI: `./a.out [gsq] [Nf] [nu0]`
- `gsq`  Wilson coupling squared (default: 12.0)
- `Nf`   number of fermion flavors, must be even (default: 2)
- `nu0`  sea mass parameter \nu_0 (default: 1.0)
- `-h`   print help and exit

Notes: ckpoint_rng files are deleted after each trajectory except every 1000th (k_ckpoint_rng=1000); ckpoint_lat files are all kept.

---

### ==meson_pq_wall_v2_claude.cu==

Partially-quenched meson propagator measurement on L=1 (N_REFINE=1, Nt=128) using wall sources with stochastic Z2 noise.
For each saved gauge configuration the program solves the overlap Dirac equation (Zolotarev order 21) with wall sources projected onto real spherical harmonics Y_{\ell m} (\ell=0..\ell_{max}, all m) and multiplied by Pauli matrices \sigma_{ab} (ab=0..3), then computes propagator inner products summed over the spatial slice.
Output is one HDF5 file per trajectory k: `data_{...}/meson_corr_v2.{k}.h5`, with datasets at path `{ell}/{em}/{ab}/{h}/{t0}/real` and `.../imag` (each a length-Nt array).
The valence mass \nu_1 (used in the Dirac operator) can differ from the sea mass \nu_0 (used to find the gauge ensemble directory), enabling partial quenching.

CLI: `./a.out [options]`
- `--gsq <gsq>`       Wilson coupling squared (default: 2.0)
- `--Nf <Nf>`         number of fermion flavors for gauge directory lookup (default: 2)
- `--nu0 <nu0>`       sea mass \nu_0 used in directory name (default: 1.0)
- `--nu1 <nu1>`       valence mass \nu_1 used in the Dirac operator (default: 1.0)
- `--nhits <nhits>`   number of stochastic hits per time slice (default: 1)
- `--dt <dt>`         spacing between wall source time slices (default: Nt/2 = 64)
- `--ellmax <ellmax>` maximum angular momentum \ell_{max} (default: 2)
- `-h, --help`        print help and exit

Notes: See tune_nu1_claude_v2.sh.

---

### meson_pq_wall_v2_L2_claude.cu

Partially-quenched meson propagator measurement for the L=2 lattice (N_REFINE=2, 42 spatial sites, Nt=128) with 4-stream parallel CG solves.
Physics and algorithm are identical to `meson_pq_wall_v2_claude.cu`, but output is plain-text files (one per trajectory, ell, em, ab, and hit) named `meson_{ell}_{em}_a{ab}_h{h}_corr_v2.{k}` with columns `t Re Im`, rather than HDF5.
Default dt is Nt=128 (a single wall source time slice t0=0 per trajectory), and k_ckpoint=1 (every trajectory is checkpointed).

CLI: `./a.out [options]`
- `--gsq <gsq>`       Wilson coupling squared (default: 2.0)
- `--Nf <Nf>`         number of fermion flavors for gauge directory lookup (default: 2)
- `--nu0 <nu0>`       sea mass \nu_0 (default: 1.0)
- `--nu1 <nu1>`       valence mass \nu_1 (default: 1.0)
- `--nhits <nhits>`   number of stochastic hits per time slice (default: 1)
- `--dt <dt>`         spacing between wall source time slices (default: Nt = 128)
- `--ellmax <ellmax>` maximum angular momentum \ell_{max} (default: 2)
- `-h, --help`        print help and exit

Notes:

---

## Jupyter Notebooks (.ipynb)

### ==meson_wall_various_ell_nu1tuning_claude_v2.ipynb==

Reads meson propagator HDF5 files produced by `meson_pq_wall_v2_claude.cu` for a list of sea/valence mass values (nu0_list), then extracts and plots the pseudoscalar (ab=3, \ell=0, m=0) correlator as a function of Euclidean time.
Compares correlator shapes for different \nu_0 values against reference exponentials exp(-\sqrt{2} t) and exp(-2t) to help locate the critical \nu_0 at which the pion mass vanishes.
Output: log-scale matplotlib plots; no output files are saved.

Notes:

---

### meson_wall_various_ell_nu1tuning.ipynb

Reads legacy plain-text meson propagator files (`meson_{ell}_{em}_a{ab}_h{h}_corr_v2.{k}`) for the (\ell=0, m=0) channel while scanning the valence mass \nu_1, with sea mass \nu_0 fixed.
Computes jackknife estimates of the ab=1,2,3 correlators summed into pseudoscalar and axial combinations (corrs_aa and corrs_33), plots the temporal correlator on a log scale, and saves a PDF `{desc}JA_corr.pdf`.
Intended for tuning the critical \nu_1 in the quenched or partially-quenched theory.

Notes:

---

### ==glue_ylms3.ipynb==

Reads glueball correlator data from `data_{...}/F_corr.{k}` and `F.{k}` files produced by `glue2_claude.cu` (16 spherical-harmonic projected plaquette operators, \ell=0..3).
Assembles the 16x16 temporal correlator matrix C(t), symmetrizes it, and solves a Generalized Eigenvalue Problem (GEVP) C(t) v = \lambda(t) C(t_0) v to extract effective masses \Delta_{eff}(t) = -\log(\lambda)/a_t for each glueball state.
Performs jackknife error analysis with configurable binsize and fits the top eigenvalue to a constant plateau to extract the lightest glueball mass; saves GEVP spectrum plot `{desc}_gevpspectrum.pdf` and fit plot `{desc}_fit_glue.pdf`.

Notes:

---

### meson_wall_various_ell.ipynb

Reads legacy plain-text meson propagator files for a scan over different (\ell, m) angular-momentum channels and Pauli matrix insertions (ab=0..3) at fixed sea and valence mass.
For each channel, computes jackknife-averaged real parts of the propagator and extracts an effective mass by fitting the ratio C(t+1)/C(t) to a single-exponential or two-state ansatz in a configurable time window [fitm, fitM].
Produces log-scale plots of correlators for multiple (\ell, m) combinations and ratio plots; saves PDF files labeled by the `desc` string.

Notes:

---

### meson_wall_v3.ipynb

Reads meson propagator HDF5 files (from `meson_pq_wall_v2_claude.cu`) for the (\ell=0, m=0, ab=0) scalar channel; extracts the effective mass meff(t) = -log(C(t+1)/C(t))/a_t and fits it to a two-state ansatz m_0 + A exp(-\Delta m t) in a user-chosen time window.
Jackknife error analysis is performed with configurable binsize, and the fitted mass is printed in the format `m(err)`.
Plots the effective mass with error bars alongside the fit curve.

Notes:

---

### meson_wall_various_ell_nu1tuning_claude.ipynb

Reads meson propagator HDF5 files from `meson_pq_wall_v2_claude.cu` for a scan over valence mass \nu_1 (nu1_list) at fixed sea mass \nu_0, extracting the (\ell=0, m=0, ab=3) pseudoscalar channel correlator.
Compares the \nu_1-dependent propagator shapes against reference exponentials exp(-\sqrt{2} t) and exp(-2t) to identify the critical valence mass, analogous to `meson_wall_various_ell_nu1tuning_claude_v2.ipynb` but scanning \nu_1 instead of \nu_0.
Output: log-scale correlator plots for each \nu_1 in the scan list; no output files saved.

Notes:

---

## Shell Scripts (.sh)

Scripts fall into three categories: SGE/SLURM job scripts (submitted via `qsub`/`sbatch`), local wrapper scripts (run directly, optionally dispatch qsub), and utility scripts.

---

### ==tune_nu1_claude.sh==

Runs `meson_pq_wall_v2_claude.o` locally on a single GPU to measure correlators while scanning the valence mass \nu_1 for \nu_1 tuning.
Active parameters: `gsq=8.0`, `Nf=2`, `nu0=1.5`, `nhits=1`, `dt=128`, `ellmax=1`.
Currently runs `nu1=nu0` on GPU 0; the commented-out batch loops for `nu1=1.5 2.0 2.5` and `nu1=1.75 2.25 2.75` show the intended scan structure.
Output via `tee gpu0.log`.

Notes:

---

### tune_nu1_claude_v2.sh

Quenched \nu_1 tuning run using `meson_pq_wall_v2_claude.o` on GPU 1.
Active parameters: `gsq=8.0`, `Nf=0` (quenched), `nu0=1.0`, `nhits=1`, `dt=128`, `ellmax=1`, `nu1=1.25`.
Structure mirrors `tune_nu1_claude.sh`; the other nu1 values are commented out.
Output via `tee gpu0.log`.

Notes: I think this was a duplicate just to run two GPU jobs in parallel

---

### run_wrapper_nf.sh

Local wrapper that runs `hmc.o` (full HMC with dynamical fermions) directly on `CUDA_VISIBLE_DEVICES=1`.
Active loop: `gsq=2.5`, `Nf=4`, `nu0=1.0`.
The qsub dispatch lines are commented out; uncomment them to submit to the BU SCC cluster instead.

Notes:

---

### ==run_wrapper_nf_v2.sh==

Same structure as `run_wrapper_nf.sh` but with different active parameters: `gsq=8.0`, `Nf=2`, `nu0=1.5`.
Runs `hmc.o` on `CUDA_VISIBLE_DEVICES=1`.

Notes:

---

### run_wrapper_glue.sh

Local wrapper that builds `glue2.o` via `make` then runs it directly for a scan over `gsq=2.0`, `Nf=2,4,6`, `nu0=1.0`.
Sources `env.sh` from the BU SCC project directory.
The qsub dispatch lines are commented out.

Notes:

---

### ==run_glue_example.sh==

Minimal one-liner example that invokes `glue2_claude.o` directly with hardcoded parameters `gsq=8.0`, `Nf=2`, `nu0=1.0`.
Intended as a quick reference for the correct argument order when running the glueball correlator measurement by hand.
No logging or looping; output goes to stdout.

Notes:

---

### run_wrapper_meson.sh

Wrapper that submits `meson.o` to the BU SCC cluster via `qsub` for a scan: `gsq=4.0`, `Nf=2,4,6`, `nu0=0.8,1.2`.
Job names are formatted `mesonNf{Nf}gsq{gsq}{nu0}`.
Uses `run_meson.sh` as the inner job script.

Notes:

---

### run_wrapper_meson_pq_wall.sh

Local wrapper that runs `meson_pq_wall.o` (old positional-arg version, not the `_claude` build) on `CUDA_VISIBLE_DEVICES=1`.
Active hardcoded parameters: `gsq=2.`, `Nf=6`, `nu0=1.0`, `nu1=nu0`, `nhits=1`, `dt=32`, `ell=0`, `em=0`.
The outer scan loops over `gsq`, `Nf`, `nu0`, `ell`, `em` are commented out.

Notes:

---

### run_wrapper.sh

Generic wrapper for the gauge-measurement binary `gauge${1}.o`.
Submits to BU SCC via qsub as an array job in the task range `${2}-${3}`, with an optional hold dependency `${4}`.
Usage: `./run_wrapper.sh <suffix> <first> <last> [hold_jid]`.

Notes:

---

### test_gauge_wrapper.sh

Wrapper that compiles `analyze_corr.cpp` (a gauge correlator analysis tool) and submits `test_gauge.sh` to the BU SCC cluster via qsub.
Currently active: `key="gsq0.500000at0.375000nt32L4_"`, `prefix_max=16`, `Nt=32`, `at=0.375`.
Many other commented-out configurations cover `L=1,2,4` at various `gsq`, `at`, `Nt` values.

Notes:

---

### run_nf_fermilab.sh

SLURM job script for the Fermilab `lq2_gpu` partition (A100 GPU).
Runs `hmc_fermilab_L2.o 12.0 6` and logs GPU utilization every 30 seconds to `gpu_usage.csv` via `nvidia-smi` in the background.
The alternative `hmc_fermilab.o` line is commented out.

Notes:

---

### run_meson_fermilab.sh

SLURM job script for the Fermilab `lq2_gpu` partition (A100 GPU).
Runs `meson_pq_wall_v2.o` with hardcoded flags `--gsq 12.0 --Nf 6 --nu0 1.0 --nu1 1.0 --nhits 1 --dt 192 --ellmax 2`.
The `meson_pq_wall_v2_L2.o` variant line is commented out.

Notes:

---

### run_nf.sh

SGE job script (BU SCC) for the full dynamical-fermion HMC.
Expects environment variables `${app}`, `${gsq}`, `${Nf}`, `${nu0}` passed via `qsub -v`.
Time limit 12 hours; requires a GPU with compute capability >= 7.0 (`gpu_c=70`).

Notes:

---

### run_meson.sh

SGE job script (BU SCC) for a (legacy) meson measurement.
Expects `${app}`, `${gsq}`, `${Nf}`, `${nu0}`.
Time limit 1 hour; V100 GPU.

Notes:

---

### run_meson_pq_wall.sh

SGE job script (BU SCC) for the partially-quenched wall-source meson measurement.
Expects `${app}`, `${gsq}`, `${Nf}`, `${nu0}`, `${nu1}`, `${ell}`, `${em}`; hardcodes `nhits=1`, `dt=32`.
Time limit 1 hour; V100 GPU.

Notes:

---

### run_glue.sh

SGE job script (BU SCC) for glueball measurements.
Expects `${app}`, `${gsq}`, `${Nf}`, `${nu0}` via `qsub -v`.
Time limit 1 hour; 1 OMP slot.

Notes:

---

### run_prop.sh

SGE job script (BU SCC) for the propagator measurement (`prop.o`).
Hardcodes `app=prop.o` with no arguments; runs with 12 OMP threads.
Time limit 1 hour.

Notes:

---

### run.sh

Generic SGE array job script (BU SCC).
Runs `${app} ${SGE_TASK_ID}`; intended for use with `-t` array task submission.
Modules: `gcc/13.2.0`, `cuda/12.5`.

Notes:

---

### test_gauge.sh

SGE job script (BU SCC) that runs `./a.out` (a gauge correlator analysis binary) with arguments from environment variables `${key}`, `${Nt}`, `${prefix_max}`, `${at}`.
The argument format is `${dir} plaq_ss_t_ 10000 10000000 50 ${Nt} 2000 ${prefix_max} ${at} "corr_${dir}"`.

Notes:

---

### remove_script.sh

Template script for deleting checkpoint and data files (`ckpoint_lat`, `ckpoint_rng`, `F`, `F_corr`) from a specific run directory over a range of trajectory indices.
All `rm` lines are commented out; edit the active `dir` and uncomment the relevant lines before running.

Notes:

---

### write.sh

One-liner utility using Ghostscript (`gs`) to extract a single page from a PDF.
Usage: `./write.sh <basename> <page>` -- reads `<basename>_.pdf`, writes `<basename>.pdf`.

Notes:

---

### zipit.sh

Archives and removes all directories matching `gsq2.000000*` in the current directory.
For each directory: creates `arxv_<dir>.tar.gz` then deletes the original with `rm -rf`.
SGE job header present but not typically submitted via qsub.

Notes:
