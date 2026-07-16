#pragma once

// multishift_kernels_claude.h
// Block kernels for the multi-shift CG (MatPoly::solve_multishift), validated host
// reference: test_multishift_claude.cpp (C1). The per-pole solutions x_m and search
// directions p_m are stored as N*npole COLUMN-MAJOR blocks [m*N + i] (m = pole index
// in [0,npole), i = lattice index in [0,N)). The per-pole CG scalars al_m, zeta_m,
// bet_m are length-npole device arrays; each thread reads ITS pole's scalar by the
// index m = gid/N ("arrays of scalars indexed by the thread id").
//
// Must be included AFTER gpu_header.h (uses CuC, Idx, NThreadsPerBlock, and the
// double*CuC / CuC+CuC operator overloads defined there).

// x_m += al_m * p_m   for all poles m, all sites i (block over N*npole).
template<Idx N> __global__
void multishift_x_update(CuC* d_x, const double* d_alm, const CuC* d_p, const int npole){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*npole;
  if(gid < total){
    const int m = gid / N;                 // pole index (N is compile-time)
    d_x[gid] = d_x[gid] + d_alm[m]*d_p[gid];
  }
}

// p_m = zeta_m * r + bet_m * p_m   for all poles m, all sites i. The seed residual
// r has length N and is broadcast over poles (read at i = gid - m*N = gid % N).
template<Idx N> __global__
void multishift_p_update(CuC* d_p, const double* d_zeta, const CuC* d_r,
                         const double* d_betm, const int npole){
  const Idx gid = (Idx)blockIdx.x*blockDim.x + threadIdx.x;
  const Idx total = (Idx)N*npole;
  if(gid < total){
    const int m = gid / N;                 // pole index
    const Idx i = gid - (Idx)m*N;          // lattice index (gid % N)
    d_p[gid] = d_zeta[m]*d_r[i] + d_betm[m]*d_p[gid];
  }
}
