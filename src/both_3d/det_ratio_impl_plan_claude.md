# det_ratio_impl_plan_claude.md

_Written: 2026-05-22 (rewrite 3)_

## Goal

Edit `eig_wmass_claude.cu` to measure the eigenspectrum of massless $D_\text{ov}$
on a single gauge configuration specified at the command line.
Output: eigenvalues (Re, Im) to an HDF5 file.
The determinant ratio is computed separately in a Python script from the HDF5 files.
A shell script loops over checkpoint numbers and calls this program once per config.

---

## Files

| File | Action |
|------|--------|
| `eig_wmass_claude.cu` | EDIT |

No new files.

---

## CLI (positional)

```
./a.out <config_dir> <k>
```

- `config_dir` — path to the gauge config directory (e.g. `data_Nf2_gsq8.0at0.03125nt128L1`)
- `k` — checkpoint number; reads `config_dir/ckpoint_lat.k`

Default behavior when args are absent: use a trivial (hot-start random) gauge field,
so the program still compiles and runs standalone for testing.

---

## Output

Write `obs_<config_dir>/eig.<k>` as a plain text file, two rows:
- Row 1: N space-separated Re($\lambda_i$) values
- Row 2: N space-separated Im($\lambda_i$) values

Output directory `obs_<config_dir>` is created if absent.

---

## Ordered chunks

### Chunk 1 — Remove mass ifdef; add filesystem include; replace CLI

- Remove `#define NO_MASS` / `#ifdef NO_MASS ... #else ... #endif` blocks.
  Always `#include "overlap.h"` and `using Fermion = Overlap<WilsonDirac>`.
  Remove `mass` variable entirely (not used anywhere).
- Add `#include <filesystem>` (for `create_directory`). No HDF5.
- Replace `atof(argv[1/2])` parsing with:

```cpp
std::string config_dir = "";
int k_cfg = -1;
if(argc > 1) config_dir = std::string(argv[1]);
if(argc > 2) k_cfg = std::atoi(argv[2]);
```

**Files:** `eig_wmass_claude.cu`

---

### Chunk 2 — Conditional config load

After `Gauge U(base)` and the existing random init:

```cpp
if(!config_dir.empty() && k_cfg >= 0){
  U.read(config_dir + "/ckpoint_lat." + std::to_string(k_cfg));
  std::cout << "# loaded " << config_dir << "/ckpoint_lat." << k_cfg << std::endl;
} else {
  std::cout << "# trivial (random) gauge field" << std::endl;
}
```

Remove `srand`/`Rng` lines since RNG is not needed for this program.

**Files:** `eig_wmass_claude.cu`

---

### Chunk 3 — Always massless operator; build matrix; run geev

No change needed to the matrix-build or geev block — it already works.
Only change: constructor is always `Fermion Dov(DW, 21)` (massless, from Chunk 1 ifdef removal).

**Files:** `eig_wmass_claude.cu`

---

### Chunk 4 — Text file output; remove old clog eigenvalue print

After `cudaMemcpy(W, d_W, CD*n, D2H)` and `assert(info==0)`:

```cpp
if(!config_dir.empty() && k_cfg >= 0){
  std::string obs_dir = "obs_" + config_dir;
  std::filesystem::create_directory(obs_dir);
  std::string outpath = obs_dir + "/eig." + std::to_string(k_cfg);
  std::ofstream ofs(outpath);
  ofs << std::scientific << std::setprecision(15);
  // row 1: Re(\lambda_i)
  for(int i = 0; i < n; i++) ofs << real(W[i]) << (i+1<n ? " " : "\n");
  // row 2: Im(\lambda_i)
  for(int i = 0; i < n; i++) ofs << imag(W[i]) << (i+1<n ? " " : "\n");
  std::cout << "# wrote " << outpath << std::endl;
} else {
  // trivial run: print to stdout
  for(int i = 0; i < n; i++)
    std::cout << real(W[i]) << " " << imag(W[i]) << std::endl;
}
```

Remove the existing `std::clog` eigenvalue print loop (line 310).

**Files:** `eig_wmass_claude.cu`

---

## Compile-time parameters

| Parameter | Debug | Production |
|-----------|-------|------------|
| `N_REFINE` | 1 | 1 |
| `Nt` | 1 (N=24) | 128 (N=3072) |

Start with `Nt=1` for compile+run check, then switch to `Nt=128`.

---

## Resolved questions

| Q | Answer |
|---|--------|
| Config loop | shell script; program handles one config |
| Mass/ratio | external Python script from HDF5 |
| Output | `obs_<config_dir>/eig.<k>`, two rows: Re then Im, space-separated |
| Trivial case | random U, print to stdout |
