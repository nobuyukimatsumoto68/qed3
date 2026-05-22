# Implementation Plan: disc_claude.cu -- Disconnected Meson Correlator

## Physics summary

The disconnected correlator is:

$$G_\text{disc}(\Delta t) = \frac{1}{N_t} \sum_{t_0} L(t_0)\cdot L^*\!\left((t_0 + \Delta t) \bmod N_t\right)$$

where the quark loop per timeslice is:

$$L(t) = \frac{1}{N_\text{sites}}\sum_{\mathbf{x}} \mathrm{Tr}\!\left[\sigma_3\, D^{-1}(t,\mathbf{x};\,t,\mathbf{x})\right]$$

Stochastic estimation: for a diluted source $\eta$,

$$L(t) \approx \frac{1}{N_\text{sites}}\sum_{\mathbf{x}} \eta^*(t,\mathbf{x},s)\,[\sigma_3]_{ss}\,\phi(t,\mathbf{x},s)$$

where $\phi = D^{-1}\eta$, and $[\sigma_3]_{ss} = +1$ for $s=0$, $-1$ for $s=1$.

Dilution: time dilution with `t_block` blocks (default $N_t/\texttt{interval} = 8$),
plus spin dilution ($N_S = 2$ inversions per block). No spatial dilution.

Fixed quantum numbers: $\ell=0$, $m=0$, ab=3 (pseudoscalar). No looping over them.

---

## Source conventions

`t_block` = number of dilution blocks (default 8 = $N_t/\texttt{interval}$).
`interval` = $N_t /$ `t_block` = spacing between active timeslices within one block (default 16).

$t_s \in [0,\, N_t/\texttt{t\_block}) = [0,\, \texttt{interval})$ labels the residue class mod `interval`.
Active timeslices for offset $t_s$: $t_s,\; t_s+\texttt{interval},\; t_s+2\cdot\texttt{interval},\;\ldots$
Complete coverage: $\texttt{interval}\times\texttt{t\_block} = N_t$.

`summed_trace_per_timeslice` is declared as:

```cpp
std::vector<Complex> summed_trace_per_timeslice(Comp::Nt, Complex(0.0, 0.0));
```

It holds one Complex scalar per timeslice:

$$L(t) = \frac{1}{N_\text{sites}}\sum_{\mathbf{x}} \eta^*(t,\mathbf{x},\text{spin})\cdot c\cdot\phi(t,\mathbf{x},\text{spin})$$

accumulated across all $(t_s, \text{spin})$ blocks into the relevant timeslices. It is NOT a `FermionVector`
(which carries full $(N_t, N_\text{sites}, N_S)$ shape); it is reset to zero at the start of each hit $h$.

For each $t_s \in [0,\, N_t/\texttt{t\_block})$ and each $\text{spin} \in \{0, 1\}$:

1. Generate diluted Z2 source: `eta.time_spin_dilution(rng, t_s, t_block, spin)`
   (nonzero only where $(t - t_s) \bmod \texttt{interval} = 0$ and spin-component $= \text{spin}$)
2. Apply $D^\dagger$: `op_DH.from_cpu<N>(DH_eta.field, eta.field)`
3. Solve $(D^\dagger D)\,\phi = D^\dagger\eta$: `op_DHD.solve<N>(phi.field, DH_eta.field)`
4. Accumulate: `eta.accumulate_loop(summed_trace_per_timeslice, phi, t_s, t_block, spin)`
   (adds $\frac{1}{N_\text{sites}}\sum_\mathbf{x}\eta^*(t,\mathbf{x},\text{spin})\cdot c\cdot\phi(t,\mathbf{x},\text{spin})$ into `summed_trace_per_timeslice[t]`,
    where $c = +1$ for $\text{spin}=0$, $c = -1$ for $\text{spin}=1$;
    the $1/N_\text{sites}$ factor (`/ Comp::N_SITES`) is applied inside `accumulate_loop`)

Total inversions per hit: $(N_t/\texttt{t\_block})\times N_S = 16\times 2 = 32$
(for default `t_block`=8, `interval`=16, $N_t=128$).

Note: `FermionVector::operator()(t, ix, spin)` returns `Complex&` (writable). Confirmed in valence.h:95.

---

## HDF5 output structure

File: `disc_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/disc_corr.{k}.h5`

Dataset path: `{h}/real` and `{h}/imag`, each a `std::vector<double>` of length $N_t$.

One file per trajectory $k$. One dataset pair per hit $h$. No $\ell$/ab nesting
(fixed to ab=3; no $t_s$ index since the time average is done inside the program).

---

## CLI (getopt_long)

Parameters (remove dt, ellmax; add t_block):

  --gsq      Wilson coupling squared (default 8.0)
  --Nf       number of fermion flavors (default 2)
  --nu0      sea quark asymmetry parameter (default 1.0)
  --nu1      valence asymmetry parameter (default 1.0; always set equal to nu0)
  --nhits    number of stochastic hits (default 1)
  --t_block  timeslices per dilution block (default 8; interval = Nt/t_block = 16)

Short opts: g, N, n, v, H, t (+ h for help).

---

## Compile-time changes (namespace Comp)

- $N_t$: 96 -> 128
- NPARALLEL block: replace the IS_OVERLAP ifdef with the flat pattern from meson:
    constexpr int NPARALLEL_DUPDATE=1;
    constexpr int NPARALLEL=NPARALLEL_DUPDATE;
    constexpr int NSTREAMS=NPARALLEL_DUPDATE;
    constexpr int NPARALLEL_GAUGE=NPARALLEL_DUPDATE;
    constexpr int NPARALLEL_SORT=NPARALLEL_DUPDATE;
- Remove commented-out $N_t$ variants (cosmetic)

---

## dir3 / dir4 naming

Dynamical ($N_f > 0$):
  dir3 = `Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nt{Nt}L{N_REFINE}/`
  dir4 = `disc_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/`

Quenched ($N_f = 0$):
  dir3 = `gsq{gsq}at{at}nt{Nt}L{N_REFINE}/`
  dir4 = `disc_gsq{gsq}at{at}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/`

($\nu_1$ appears in dir4 only; dir3 is identical to the meson code.)

---

## Checkpoint / resume logic

kmax = 1000 (unchanged from meson file).
k_ckpoint = 10.

k_tmp scan: walk $k = 10, 20, \ldots$ checking `ckpoint_lat.{k}` exists; stop at last found.

k_init resume scan: walk $k = 10, 20, \ldots$ checking `disc_corr.{k}.h5` exists; stop at first
missing file, then back up one step.

Resume check inside the $k$ loop: if `disc_corr.{k}.h5` exists, open ReadOnly and check
  `h5check.exist( std::to_string(nhits-1) + "/real" )`
If true, skip (continue).

---

## Ordered implementation chunks

### Chunk 1 -- Headers and namespace Comp

**Files:** `disc_claude.cu`

- Add `#include <highfive/H5File.hpp>` and `#include <getopt.h>`
- Remove `#include "../../geometry/geodesic.h"` (line 99)
- Change `const std::string dir = "../../dats/";` -> `"../../geometry/data/"`
- Fix $N_t$: 96 -> 128
- Replace NPARALLEL block (IS_OVERLAP ifdef, lines 46-58) with flat NPARALLEL_DUPDATE=1 pattern
- Remove commented-out $N_t$ variants

### Chunk 2 -- PrintHelp() and ParseArgs()

**Files:** `disc_claude.cu`

Insert before main(). Pattern copied from meson_pq_wall_v2_claude.cu, with:
- Parameters: gsq, Nf, nu0, nu1, nhits, t_block
- t_block: short opt 't', long opt "t_block", default 8
- Remove: dt, ellmax

### Chunk 3 -- main() preamble

**Files:** `disc_claude.cu`

Replace positional-arg parsing and igam block with:
  `ParseArgs(argc, argv, gsq, Nf, nu0, nu1, nhits, t_block);`
Print summary to stdout.

### Chunk 4 -- dir3 / dir4

**Files:** `disc_claude.cu`

Copy the quenched/dynamical branching pattern from meson; change dir4 prefix
from `data_Nf...` to `disc_Nf...` (or `disc_gsq...` for $N_f=0$).
Add $\nu_1$ and `tb{t_block}` to the dir4 name. `std::filesystem::create_directory(dir4)`.

### Chunk 5 -- Dirac operator setup

**Files:** `disc_claude.cu`

- Remove IS_DUAL branch entirely (keep S2Simp / DiracS2Simp)
- Remove IS_DAGGER define usage
- Remove IS_OVERLAP #else branch (DiracPf path, lines 277-285)
- Change WilsonDirac DW asymmetry argument: $\nu_0$ -> $\nu_1$ (use valence asymmetry parameter; $\nu_1 = \nu_0$ in practice)
- Keep `Fermion D(DW, 21)`
- Remove op_D and op_DDH; keep op_DH and op_DHD only

### Chunk 6 -- FermionVector declarations

**Files:** `disc_claude.cu`

Remove: src0, src1, DHsrc0, DHsrc1, Dinv0, Dinv1, gauge (GAUGE_TRSF block).
Declare:
  `FermionVector eta;`
  `FermionVector DH_eta;`
  `FermionVector phi;`
  `std::vector<Complex> summed_trace_per_timeslice(Comp::Nt, Complex(0.0, 0.0));`
  `std::vector<double> G_real(Comp::Nt, 0.0), G_imag(Comp::Nt, 0.0);`

### Chunk 7 -- Checkpoint scan and k_init resume scan

**Files:** `disc_claude.cu`

Copy from meson_pq_wall_v2_claude.cu (lines 400-423), with:
- filename `disc_corr.` instead of `meson_corr_v2.`
- resume sentinel: `std::to_string(nhits-1) + "/real"`

### Chunk 8 -- Trajectory loop (outer structure)

**Files:** `disc_claude.cu`

```cpp
for(int k=k_init; k<=k_tmp; k+=k_ckpoint){
    // --- resume check ---
    const std::string h5path = dir4+"disc_corr."+std::to_string(k)+".h5";
    if(std::filesystem::exists(h5path)){
      try {
        HighFive::File h5check(h5path, HighFive::File::ReadOnly);
        const std::string last_ds = std::to_string(nhits-1)+"/real";
        if(h5check.exist(last_ds)){
          std::cout << "# skipping k = " << k << " (complete)" << std::endl;
          continue;
        }
      } catch(...) {}
    }
    // --- read gauge config and update Dirac op ---
    std::cout << "# k = " << k << std::endl;
    U.read( dir3+"ckpoint_lat."+std::to_string(k) );
    D.update( U );
    // --- open HDF5 for this trajectory ---
    HighFive::File h5file(h5path,
        HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    // --- stochastic hits ---
    for(int h=0; h<nhits; h++){
        // ... Chunk 9: fill summed_trace_per_timeslice ...
        // ... Chunk 10: compute G_real, G_imag ...
        // ... Chunk 11: write HDF5 ...
    }
} // end for k
```

Note: Truncate means a partially-written file (interrupted run) is overwritten from scratch
on restart. The k_init resume scan skips only fully-complete files.

### Chunk 9 -- Stochastic loop (inside for h)

**Files:** `valence.h` (add `time_spin_dilution`, `accumulate_loop`), `disc_claude.cu`

```cpp
std::fill(summed_trace_per_timeslice.begin(), summed_trace_per_timeslice.end(), Complex(0.0, 0.0));

const int interval = Comp::Nt / t_block;
for(int t_s=0; t_s<interval; t_s++){
    for(int spin=0; spin<NS; spin++){
        eta.time_spin_dilution(rng, t_s, t_block, spin);
        op_DH.from_cpu<N>(DH_eta.field, eta.field);
        op_DHD.solve<N>(phi.field, DH_eta.field);
        eta.accumulate_loop(summed_trace_per_timeslice, phi, t_s, t_block, spin);
    }
}
```

### Chunk 10 -- Correlator computation (inside for h, after Chunk 9)

**Files:** `disc_claude.cu` (add `compute_disc_corr` free function before `main`)

Computes:

$$G_\text{disc}[\Delta t] = \frac{1}{N_t}\sum_{t_0} L[t_0]\cdot L^*\!\left[(t_0+\Delta t)\bmod N_t\right]$$

```cpp
std::fill(G_real.begin(), G_real.end(), 0.0);
std::fill(G_imag.begin(), G_imag.end(), 0.0);

for(int t0=0; t0<Comp::Nt; t0++){
    for(int dt=0; dt<Comp::Nt; dt++){
        const Complex contrib = summed_trace_per_timeslice[t0]
            * std::conj(summed_trace_per_timeslice[(t0+dt) % Comp::Nt]);
        G_real[dt] += contrib.real() / Comp::Nt;
        G_imag[dt] += contrib.imag() / Comp::Nt;
    }
}
```

### Chunk 11 -- HDF5 write (inside for h, after Chunk 10)

**Files:** `disc_claude.cu`

```cpp
const std::string ds_prefix = std::to_string(h) + "/";
h5file.createDataSet(ds_prefix + "real", G_real);
h5file.createDataSet(ds_prefix + "imag", G_imag);
```

### Chunk 12 -- Remove old free-field block

**Files:** `disc_claude.cu`

Remove entirely: the igam loop (lines 360-382) and all old FermionVector
uses (src0, src1, Dinv0, Dinv1, etc.), and the commented-out gauge-config loop (lines 387-447).
Keep: `for(int i=0; i<Comp::NSTREAMS; i++) d_MemorySets[i].deallocate();`

---

## Open questions (resolved)

1. ~~**DW asymmetry argument (Chunk 5):**~~ **RESOLVED.** $\nu_0$ is the asymmetry parameter
   (speed of light), not the sea quark mass. $\nu_1$ is the valence asymmetry parameter;
   $\nu_1 = \nu_0$ always unless otherwise noted. Use `WilsonDirac DW(..., nu1)`.

2. ~~**Chunk 9 accumulation:**~~ **RESOLVED.** Accumulate $L(t)$ only at active timeslices
   ($(t - t_s)\bmod\texttt{interval} = 0$). Standard interlaced time-dilution estimator.

3. ~~**Additional routines needed?**~~ **RESOLVED.** CPU loops sufficient; no GPU kernel needed.
   `time_spin_dilution` and `accumulate_loop` are CPU methods in `valence.h`.

4. ~~**Output directory name:**~~ **RESOLVED.** Append `tb{t_block}` to dir4:
   `disc_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/`
   (quenched: `disc_gsq{gsq}at{at}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/`)

---

## Validation checks

### C1 -- Correlator symmetry (coding check, exact)

From the $t_0$ sum: $G(\Delta t) = G^*(N_t - \Delta t)$.
Therefore:

- $G_\text{real}(\Delta t) = G_\text{real}(N_t - \Delta t)$ (even)
- $G_\text{imag}(\Delta t) = -G_\text{imag}(N_t - \Delta t)$ (odd)

These hold exactly for any $L(t)$, regardless of statistics. Any violation is a coding bug.

Corollary (from oddness at fixed points): $G_\text{imag}(0) = 0$ and $G_\text{imag}(N_t/2) = 0$ exactly.

### C2 -- Positivity at $\Delta t = 0$ (exact)

$$G_\text{real}(0) = \frac{1}{N_t}\sum_{t_0} |L(t_0)|^2 \geq 0$$

always. Should be positive for any nontrivial gauge config.

### C3 -- Nhits convergence (stochastic check)

For a fixed gauge config, run with `nhits` = 1, 2, 4, 8.
The sample mean of $G_\text{real}(\Delta t)$ should be stable; the standard deviation across hits
should scale as $1/\sqrt{\text{nhits}}$.

### C4 -- Free-field time-translation invariance

Initialize $U = 1$ (all gauge links zero) by skipping the file read and calling `D.update(U)` on a
fresh `Gauge` object. Run with many hits. Then $|L(t)|$ should be independent of $t$ to within
stochastic noise, since the free operator is time-translation invariant.

### C5 -- Large-$\Delta t$ exponential decay

For a confining or massive phase, $G_\text{disc}(\Delta t) \sim A\,e^{-m_\eta \Delta t}$ at large $\Delta t$.
The extracted mass $m_\eta$ should be consistent with the eta-prime mass measured by other methods
(connected + disconnected). This is a physics cross-check, not a coding check.

### C6 -- Reality of $L(t)$ in the free field (free-field only)

In the free-field limit ($U=1$), the 3D two-component overlap operator admits $\sigma_3$ as a
Hermiticity matrix, so $\sigma_3 D \sigma_3 = D^\dagger$ holds and the argument from the
former C3 applies: $L(t)$ is real. Therefore $G_\text{imag}(\Delta t) = 0$ exactly.

Check: run C4 (free-field, $U=1$) and verify $|G_\text{imag}(\Delta t)|$ is consistent with zero
within stochastic noise. This is stronger than C1 and provides an additional coding/physics check
specific to the free-field configuration.
