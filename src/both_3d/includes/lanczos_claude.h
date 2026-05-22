#pragma once
// lanczos_claude.h
// Implicit Restarted Lanczos (IRL) eigensolver with Chebyshev polynomial
// acceleration for the operator A = (D_ov + m)^\dagger (D_ov + m).
//
// Parameters:
//   Nk  = 50   converged eigenpairs sought
//   Nm  = 150  total Krylov dimension
//   \alpha  = upper edge of target window (smaller value)
//   \beta   = spectral upper bound (larger value)
//   cheb_degree = degree of Chebyshev polynomial (default 12)
//
// Chebyshev filter (Sec 3.2 of lanczos.pdf):
//   q(\mu) = (2\mu - (\alpha^2 + \beta^2)) / (\alpha^2 - \beta^2)
//   Applied via T_d(q(A)) using the three-term recurrence.
//
// Synchronization rule: after every Taxpy or Taxpy_gen kernel launch,
// call CUDA_CHECK(cudaDeviceSynchronize()) unconditionally.

#include <vector>
#include <cmath>
#include <Eigen/Dense>

// ============================================================
// Chunk 1: ChebyshevFilter
// ============================================================

// ChebyshevFilter applies T_d(q(A)) v via the three-term recurrence
//   p_0 = v
//   p_1 = q(A) v
//   p_{j+1} = 2 q(A) p_j - p_{j-1}
// where q(\mu) = (2\mu - (\alpha^2 + \beta^2)) / (\alpha^2 - \beta^2).
// \alpha < \beta: alpha = upper edge of target window, beta = spectral upper bound.
struct ChebyshevFilter {
  const LinOp& _Op;   // HermOp: A = (D+m)^\dagger(D+m)
  const int    _deg;  // polynomial degree d
  const double _alpha; // upper edge of target window (smaller)
  const double _beta;  // spectral upper bound (larger)

  // Scratch device buffers (owned by this object).
  // Declared non-const so we can initialize them in constructor.
  CuC* _d_v0;
  CuC* _d_v1;
  CuC* _d_tmp;

  ChebyshevFilter(const LinOp& Op,
                  const int deg,
                  const double alpha,
                  const double beta)
    : _Op(Op)
    , _deg(deg)
    , _alpha(alpha)
    , _beta(beta)
  {
    CUDA_CHECK(cudaMalloc(&_d_v0,  Comp::N * CD));
    CUDA_CHECK(cudaMalloc(&_d_v1,  Comp::N * CD));
    CUDA_CHECK(cudaMalloc(&_d_tmp, Comp::N * CD));
  }

  ~ChebyshevFilter()
  {
    CUDA_CHECK(cudaFree(_d_v0));
    CUDA_CHECK(cudaFree(_d_v1));
    CUDA_CHECK(cudaFree(_d_tmp));
  }

  // operator()(d_out, d_in): apply T_d(q(A)) to d_in, write result to d_out.
  // Uses local pointer variables to avoid mutating const members.
  void operator()(CuC* d_out, const CuC* d_in) const
  {
    constexpr Idx N = Comp::N;

    // q(A) v = (2 A v - (\alpha^2 + \beta^2) v) / (\alpha^2 - \beta^2)
    // Precompute scalar coefficients.
    const double a2 = _alpha * _alpha;
    const double b2 = _beta  * _beta;
    const double denom = a2 - b2;           // negative: \alpha < \beta
    const double shift = -(a2 + b2) / denom; // = (\alpha^2 + \beta^2) / (\beta^2 - \alpha^2)
    const double scale =  2.0 / denom;       // = -2 / (\beta^2 - \alpha^2)

    // Use local aliases so we can swap without const issues.
    CuC* p0 = _d_v0;
    CuC* p1 = _d_v1;
    CuC* tmp = _d_tmp;

    // p0 = v (copy d_in into p0)
    CUDA_CHECK(cudaMemcpy(p0, d_in, N * CD, D2D));

    if (_deg == 0) {
      CUDA_CHECK(cudaMemcpy(d_out, p0, N * CD, D2D));
      return;
    }

    // p1 = q(A) p0 = scale * A*p0 + shift * p0
    _Op(tmp, p0);   // tmp = A p0; _Op calls cudaDeviceSynchronize inside
    // p1 = scale * tmp  (zero p1 first)
    CUDA_CHECK(cudaMemset(p1, 0, N * CD));
    Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
      p1, scale, tmp, p1);
    CUDA_CHECK(cudaDeviceSynchronize());
    // p1 += shift * p0
    Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
      p1, shift, p0, p1);
    CUDA_CHECK(cudaDeviceSynchronize());

    if (_deg == 1) {
      CUDA_CHECK(cudaMemcpy(d_out, p1, N * CD, D2D));
      return;
    }

    // Three-term recurrence: p_{j+1} = 2 q(A) p_j - p_{j-1}
    // We need a third buffer; use d_out as the output/next buffer.
    CuC* pnew = d_out;

    for (int j = 1; j < _deg; ++j) {
      // tmp2 = q(A) p1 = scale * A p1 + shift * p1
      _Op(tmp, p1);  // tmp = A p1
      // pnew = 2*scale * tmp + 2*shift * p1 - p0
      // Step: pnew = 2*scale * tmp
      CUDA_CHECK(cudaMemset(pnew, 0, N * CD));
      const double two_scale = 2.0 * scale;
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        pnew, two_scale, tmp, pnew);
      CUDA_CHECK(cudaDeviceSynchronize());
      // pnew += 2*shift * p1
      const double two_shift = 2.0 * shift;
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        pnew, two_shift, p1, pnew);
      CUDA_CHECK(cudaDeviceSynchronize());
      // pnew += -1.0 * p0
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        pnew, -1.0, p0, pnew);
      CUDA_CHECK(cudaDeviceSynchronize());

      // Rotate: p0 <- p1, p1 <- pnew
      // For the next iteration pnew will be reused (it points to d_out
      // which is fine as long as j < _deg - 1; on the last pass the
      // result stays in d_out).
      // Swap p0 and p1, then set pnew to the old p0.
      CuC* ptmp = p0;
      p0   = p1;
      p1   = pnew;
      pnew = ptmp;
    }

    // After the loop p1 holds the result T_d(q(A)) v.
    // If p1 != d_out, copy.
    if (p1 != d_out) {
      CUDA_CHECK(cudaMemcpy(d_out, p1, N * CD, D2D));
    }
  }
};


// ============================================================
// Chunk 2: IRLState
// ============================================================

// IRLState holds all Krylov basis vectors and tridiagonal data
// for the Implicitly Restarted Lanczos algorithm.
struct IRLState {
  const int Nk;  // number of converged eigenpairs sought
  const int Nm;  // total Krylov dimension

  // Device Krylov basis vectors: Nm vectors each of length Comp::N.
  std::vector<CuC*> d_evec;

  // Device residual vector f (length Comp::N).
  CuC* d_f;

  // Scratch buffer for basisRotate (length Comp::N).
  CuC* d_tmp_rot;

  // Tridiagonal data (host).
  std::vector<double> lmd;  // diagonal elements, length Nm
  std::vector<double> lme;  // sub-diagonal elements, length Nm

  // Converged eigenvalues.
  std::vector<double> eval;

  // Rotation matrix Qt, Nm x Nm.
  Eigen::MatrixXd Qt;

  IRLState(const int Nk_, const int Nm_)
    : Nk(Nk_)
    , Nm(Nm_)
    , d_evec(Nm_, nullptr)
    , lmd(Nm_, 0.0)
    , lme(Nm_, 0.0)
    , Qt(Nm_, Nm_)
  {
    for (int i = 0; i < Nm_; ++i) {
      CUDA_CHECK(cudaMalloc(&d_evec[i], Comp::N * CD));
      CUDA_CHECK(cudaMemset(d_evec[i], 0, Comp::N * CD));
    }
    CUDA_CHECK(cudaMalloc(&d_f,       Comp::N * CD));
    CUDA_CHECK(cudaMalloc(&d_tmp_rot, Comp::N * CD));
    CUDA_CHECK(cudaMemset(d_f,       0, Comp::N * CD));
    CUDA_CHECK(cudaMemset(d_tmp_rot, 0, Comp::N * CD));
  }

  ~IRLState()
  {
    for (int i = 0; i < Nm; ++i) {
      CUDA_CHECK(cudaFree(d_evec[i]));
    }
    CUDA_CHECK(cudaFree(d_f));
    CUDA_CHECK(cudaFree(d_tmp_rot));
  }
};


// ============================================================
// Chunk 3: lanczos_step
// ============================================================

// lanczos_step performs one step of the Lanczos iteration:
//   w = PolyOp * evec[k]
//   if k > 0: w -= lme[k-1] * evec[k-1]
//   lmd[k] = Re(<evec[k], w>)
//   w -= lmd[k] * evec[k]
//   lme[k] = ||w||
//   evec[k+1] = w / lme[k]
// Re-orthogonalization against all previous vectors every orth_period steps.
static void lanczos_step(const ChebyshevFilter& PolyOp,
                         IRLState& s,
                         const int k,
                         cublasHandle_t handle)
{
  constexpr Idx N = Comp::N;
  constexpr int orth_period = 4;

  // w = PolyOp(evec[k]) stored in s.d_f
  PolyOp(s.d_f, s.d_evec[k]);

  // w -= lme[k-1] * evec[k-1]
  if (k > 0) {
    const double neg_beta = -s.lme[k - 1];
    Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
      s.d_f, neg_beta, s.d_evec[k - 1], s.d_f);
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // lmd[k] = Re(<evec[k], w>) via cublasZdotc
  CuC dot_val;
  CUBLAS_CHECK(cublasZdotc(handle, N,
                            s.d_evec[k], 1,
                            s.d_f, 1,
                            &dot_val));
  s.lmd[k] = cuCreal(dot_val);

  // w -= lmd[k] * evec[k]
  const double neg_alpha = -s.lmd[k];
  Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
    s.d_f, neg_alpha, s.d_evec[k], s.d_f);
  CUDA_CHECK(cudaDeviceSynchronize());

  // Re-orthogonalize against all previous vectors every orth_period steps.
  if ((k > 0) && ((k % orth_period) == 0)) {
    for (int j = 0; j <= k; ++j) {
      CuC c;
      CUBLAS_CHECK(cublasZdotc(handle, N,
                                s.d_evec[j], 1,
                                s.d_f, 1,
                                &c));
      // f -= c * evec[j]  (subtract projection)
      // Taxpy_gen: d_res = d_a * d_x + d_y with T2=CuC
      Taxpy<CuC, N><<<NBlocks, NThreadsPerBlock>>>(
        s.d_f, -c, s.d_evec[j], s.d_f);
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }

  // lme[k] = ||w|| = sqrt(dot(w, w))
  CuC norm_sq;
  CUBLAS_CHECK(cublasZdotc(handle, N,
                            s.d_f, 1,
                            s.d_f, 1,
                            &norm_sq));
  const double beta = std::sqrt(cuCreal(norm_sq));
  s.lme[k] = beta;

  // evec[k+1] = w / beta   (only if k+1 < Nm)
  if (k + 1 < s.Nm) {
    if (beta > 1.0e-20) {
      const double inv_beta = 1.0 / beta;
      CUBLAS_CHECK(cublasZdscal(handle, N,
                                 &inv_beta,
                                 s.d_f, 1));
    }
    CUDA_CHECK(cudaMemcpy(s.d_evec[k + 1], s.d_f, N * CD, D2D));
  }
}


// ============================================================
// Chunk 4: diagonalize_Eigen
// ============================================================

// diagonalize_Eigen: CPU tridiagonal eigensolver using Eigen.
// On return, lmd[0..Nk-1] holds eigenvalues in DESCENDING order,
// and Qt rows 0..Nk-1 hold the corresponding eigenvectors.
// Direct port from Grid IRL (lines 501-522).
static void diagonalize_Eigen(std::vector<double>& lmd,
                               std::vector<double>& lme,
                               int Nk, int Nm,
                               Eigen::MatrixXd& Qt)
{
  Eigen::MatrixXd TriDiag = Eigen::MatrixXd::Zero(Nk, Nk);
  for (int i = 0; i < Nk; ++i)     TriDiag(i, i)     = lmd[i];
  for (int i = 0; i < Nk - 1; ++i) TriDiag(i, i + 1) = lme[i];
  for (int i = 0; i < Nk - 1; ++i) TriDiag(i + 1, i) = lme[i];

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(TriDiag);

  // Store in descending order (eigenvalues() returns ascending).
  for (int i = 0; i < Nk; ++i) {
    lmd[Nk - 1 - i] = eig.eigenvalues()(i);
  }
  for (int i = 0; i < Nk; ++i) {
    for (int j = 0; j < Nk; ++j) {
      Qt(Nk - 1 - i, j) = eig.eigenvectors()(j, i);
    }
  }
}


// ============================================================
// Chunk 5: QR_decomp
// ============================================================

// QR_decomp: one step of the implicitly shifted QR algorithm.
// Applies a sequence of Givens rotations to the tridiagonal matrix
// (lmd, lme) with shift Dsh, accumulating the rotation in Qt.
// Direct port from Grid IRL (lines 527-585).
static void QR_decomp(std::vector<double>& lmd,
                       std::vector<double>& lme,
                       int Nk, int Nm,
                       Eigen::MatrixXd& Qt,
                       double Dsh,
                       int kmin, int kmax)
{
  int k = kmin - 1;

  double Fden = 1.0 / std::hypot(lmd[k] - Dsh, lme[k]);
  double c =  (lmd[k] - Dsh) * Fden;
  double s = -lme[k] * Fden;

  double tmpa1 = lmd[k];
  double tmpa2 = lmd[k + 1];
  double tmpb  = lme[k];

  lmd[k]     = c*c*tmpa1 + s*s*tmpa2 - 2.0*c*s*tmpb;
  lmd[k + 1] = s*s*tmpa1 + c*c*tmpa2 + 2.0*c*s*tmpb;
  lme[k]     = c*s*(tmpa1 - tmpa2) + (c*c - s*s)*tmpb;
  double x   = -s * lme[k + 1];
  lme[k + 1] =  c * lme[k + 1];

  for (int i = 0; i < Nk; ++i) {
    double Qtmp1 = Qt(k, i);
    double Qtmp2 = Qt(k + 1, i);
    Qt(k,     i) =  c*Qtmp1 - s*Qtmp2;
    Qt(k + 1, i) =  s*Qtmp1 + c*Qtmp2;
  }

  // Givens transformations for the remaining rows.
  for (int kk = kmin; kk < kmax - 1; ++kk) {
    Fden = 1.0 / std::hypot(x, lme[kk - 1]);
    c =  lme[kk - 1] * Fden;
    s = -x * Fden;

    tmpa1 = lmd[kk];
    tmpa2 = lmd[kk + 1];
    tmpb  = lme[kk];

    lmd[kk]     = c*c*tmpa1 + s*s*tmpa2 - 2.0*c*s*tmpb;
    lmd[kk + 1] = s*s*tmpa1 + c*c*tmpa2 + 2.0*c*s*tmpb;
    lme[kk]     = c*s*(tmpa1 - tmpa2) + (c*c - s*s)*tmpb;
    lme[kk - 1] = c*lme[kk - 1] - s*x;

    if (kk != kmax - 2) {
      x         = -s * lme[kk + 1];
      lme[kk + 1] =  c * lme[kk + 1];
    }

    for (int i = 0; i < Nk; ++i) {
      double Qtmp1 = Qt(kk,     i);
      double Qtmp2 = Qt(kk + 1, i);
      Qt(kk,     i) =  c*Qtmp1 - s*Qtmp2;
      Qt(kk + 1, i) =  s*Qtmp1 + c*Qtmp2;
    }
  }
}


// ============================================================
// Chunk 6: basisRotate and basisRotateJ
// ============================================================

// basisRotate: rotate d_evec[k_lo..k_hi) by Qt.
// For each output index j in [k_lo, k_hi):
//   new_evec[j] = sum_{i=0}^{Nm-1} Qt(i,j) * d_evec[i]
// All outputs are accumulated into temporary device buffers first,
// then copied back, to avoid reading overwritten source vectors.
// d_tmp is a single scratch buffer used for accumulation (one j at a time).
static void basisRotate(std::vector<CuC*>& d_evec,
                         Eigen::MatrixXd& Qt,
                         int k_lo, int k_hi, int Nm,
                         CuC* d_tmp,
                         cublasHandle_t handle)
{
  constexpr Idx N = Comp::N;
  (void)handle;

  const int n_out = k_hi - k_lo;
  if (n_out <= 0) return;

  // Allocate temporary device buffers for all output vectors.
  std::vector<CuC*> d_out(n_out, nullptr);
  for (int jj = 0; jj < n_out; ++jj) {
    CUDA_CHECK(cudaMalloc(&d_out[jj], N * CD));
  }

  // Compute each output vector into d_out[jj].
  for (int jj = 0; jj < n_out; ++jj) {
    const int j = k_lo + jj;
    CUDA_CHECK(cudaMemset(d_out[jj], 0, N * CD));

    for (int i = 0; i < Nm; ++i) {
      const double coeff = Qt(j, i);
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        d_out[jj], coeff, d_evec[i], d_out[jj]);
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }

  // Copy results back into d_evec.
  for (int jj = 0; jj < n_out; ++jj) {
    const int j = k_lo + jj;
    CUDA_CHECK(cudaMemcpy(d_evec[j], d_out[jj], N * CD, D2D));
    CUDA_CHECK(cudaFree(d_out[jj]));
  }
}


// basisRotateJ: rotate a single vector j.
//   d_out = sum_{i=i_lo}^{i_hi-1} Qt(i, j) * d_evec[i]
// (Grid calls this basisRotateJ for convergence testing.)
static void basisRotateJ(CuC* d_out,
                          std::vector<CuC*>& d_evec,
                          Eigen::MatrixXd& Qt,
                          int j, int i_lo, int i_hi, int Nm,
                          CuC* d_tmp)
{
  constexpr Idx N = Comp::N;
  (void)d_tmp; // not needed; d_out serves as accumulator

  CUDA_CHECK(cudaMemset(d_out, 0, N * CD));

  for (int i = i_lo; i < i_hi; ++i) {
    const double coeff = Qt(j, i);
    Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
      d_out, coeff, d_evec[i], d_out);
    CUDA_CHECK(cudaDeviceSynchronize());
  }
}


// ============================================================
// Chunk 7: calc -- main IRL driver
// ============================================================

// calc: run the Implicitly Restarted Lanczos algorithm.
//
// Arguments:
//   PolyOp   -- Chebyshev-filtered operator applied during Lanczos steps
//   HermOp   -- original Hermitian operator A = (D+m)^\dagger(D+m)
//   s        -- IRLState (Nk, Nm already set)
//   d_src    -- initial device vector (not necessarily normalized)
//   eresid   -- convergence tolerance on residual
//   MaxIter  -- maximum number of restart iterations
//   handle   -- cuBLAS handle
static void calc(const ChebyshevFilter& PolyOp,
                 const LinOp& HermOp,
                 IRLState& s,
                 const CuC* d_src,
                 const double eresid,
                 const int MaxIter,
                 cublasHandle_t handle)
{
  constexpr Idx N = Comp::N;
  const int Nk = s.Nk;
  const int Nm = s.Nm;
  const int k1 = 1;
  const int k2 = Nk;

  // ----------------------------------------------------------
  // Power iteration: estimate largest eigenvalue of HermOp.
  // ----------------------------------------------------------
  double evalMaxApprox = 0.0;
  {
    CuC* d_src_n;
    CuC* d_tmp_pi;
    CUDA_CHECK(cudaMalloc(&d_src_n,  N * CD));
    CUDA_CHECK(cudaMalloc(&d_tmp_pi, N * CD));
    CUDA_CHECK(cudaMemcpy(d_src_n, d_src, N * CD, D2D));

    for (int i = 0; i < 50; ++i) {
      // Normalize d_src_n.
      CuC nn_c;
      CUBLAS_CHECK(cublasZdotc(handle, N,
                                d_src_n, 1,
                                d_src_n, 1,
                                &nn_c));
      const double nn = std::sqrt(cuCreal(nn_c));
      if (nn > 1.0e-20) {
        const double inv_nn = 1.0 / nn;
        CUBLAS_CHECK(cublasZdscal(handle, N, &inv_nn, d_src_n, 1));
      }

      // Apply HermOp.
      HermOp(d_tmp_pi, d_src_n);

      // Estimate: na = sqrt(||H v||^2 / ||v||^2) = ||H v|| (since v is normalized).
      CuC vnum_c;
      CUBLAS_CHECK(cublasZdotc(handle, N,
                                d_tmp_pi, 1,
                                d_tmp_pi, 1,
                                &vnum_c));
      CuC vden_c;
      CUBLAS_CHECK(cublasZdotc(handle, N,
                                d_src_n, 1,
                                d_src_n, 1,
                                &vden_c));
      const double na = std::sqrt(cuCreal(vnum_c) / cuCreal(vden_c));
      if (std::fabs(evalMaxApprox / (na + 1.0e-300) - 1.0) < 0.0001)
        break;
      evalMaxApprox = na;

      // d_src_n = d_tmp_pi (for next iteration).
      CUDA_CHECK(cudaMemcpy(d_src_n, d_tmp_pi, N * CD, D2D));
    }

    std::clog << "# evalMaxApprox = " << evalMaxApprox << std::endl;

    CUDA_CHECK(cudaFree(d_src_n));
    CUDA_CHECK(cudaFree(d_tmp_pi));
  }

  // ----------------------------------------------------------
  // Step 1: Normalize d_src into s.d_evec[0].
  // ----------------------------------------------------------
  CUDA_CHECK(cudaMemcpy(s.d_evec[0], d_src, N * CD, D2D));
  {
    CuC nn_c;
    CUBLAS_CHECK(cublasZdotc(handle, N,
                              s.d_evec[0], 1,
                              s.d_evec[0], 1,
                              &nn_c));
    const double inv_nn = 1.0 / std::sqrt(cuCreal(nn_c));
    CUBLAS_CHECK(cublasZdscal(handle, N, &inv_nn, s.d_evec[0], 1));
  }

  // ----------------------------------------------------------
  // Step 2: Run Nk initial Lanczos steps with PolyOp.
  // ----------------------------------------------------------
  for (int k = 0; k < Nk; ++k) {
    lanczos_step(PolyOp, s, k, handle);
  }
  std::clog << "# Initial " << Nk << " Lanczos steps done." << std::endl;

  // ----------------------------------------------------------
  // Step 3: Restart loop.
  // ----------------------------------------------------------
  // Working copies for diagonalize_Eigen (avoids overwriting s.lmd/s.lme).
  std::vector<double> eval2(Nm, 0.0);
  std::vector<double> lme2(Nm, 0.0);
  std::vector<double> eval2_copy(Nm, 0.0);

  // Number of shifts to apply = Nm - Nk.
  const int Np = Nm - Nk;

  bool converged_flag = false;
  int iter = 0;

  for (iter = 0; iter < MaxIter; ++iter) {
    std::clog << "# Restart iteration " << iter << std::endl;

    // a) Extend from Nk to Nm steps.
    for (int k = Nk; k < Nm; ++k) {
      lanczos_step(PolyOp, s, k, handle);
    }
    // Set d_f = lme[Nm-1] * evec[Nm-1]  (Grid: f *= lme[Nm-1]).
    {
      const double sc = s.lme[Nm - 1];
      CUDA_CHECK(cudaMemcpy(s.d_f, s.d_evec[Nm - 1], N * CD, D2D));
      CUBLAS_CHECK(cublasZdscal(handle, N, &sc, s.d_f, 1));
    }

    // b) Diagonalize the Nm x Nm tridiagonal matrix.
    for (int k = 0; k < Nm; ++k) {
      eval2[k] = s.lmd[k + k1 - 1];
      lme2[k]  = s.lme[k + k1 - 1];
    }
    s.Qt = Eigen::MatrixXd::Identity(Nm, Nm);
    diagonalize_Eigen(eval2, lme2, Nm, Nm, s.Qt);
    std::clog << "# Diagonalized." << std::endl;

    // Sort descending; largest Np eigenvalues are the shifts.
    eval2_copy = eval2;
    std::partial_sort(eval2.begin(), eval2.begin() + Nm, eval2.end(),
                      std::greater<double>());

    // c) Reset Qt and apply Np QR shifts (indices k2..Nm-1 in sorted eval2).
    s.Qt = Eigen::MatrixXd::Identity(Nm, Nm);
    for (int ip = k2; ip < Nm; ++ip) {
      QR_decomp(s.lmd, s.lme, Nm, Nm, s.Qt, eval2[ip], k1, Nm);
    }

    // d) Rotate the basis vectors: indices 0..k2 (Grid: k1-1=0, k2+1=k2+1).
    basisRotate(s.d_evec, s.Qt, 0, k2 + 1, Nm, s.d_tmp_rot, handle);

    // e) Update residual vector and lme[k2-1].
    //    f_new = Qt(k2-1, Nm-1) * d_f + lme[k2-1] * evec[k2]
    {
      const double qt_coeff  = s.Qt(k2 - 1, Nm - 1);
      const double lme_coeff = s.lme[k2 - 1];

      // d_tmp_rot = qt_coeff * d_f
      CUDA_CHECK(cudaMemset(s.d_tmp_rot, 0, N * CD));
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        s.d_tmp_rot, qt_coeff, s.d_f, s.d_tmp_rot);
      CUDA_CHECK(cudaDeviceSynchronize());
      // d_tmp_rot += lme_coeff * evec[k2]
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        s.d_tmp_rot, lme_coeff, s.d_evec[k2], s.d_tmp_rot);
      CUDA_CHECK(cudaDeviceSynchronize());

      // beta_k = ||d_tmp_rot||
      CuC norm_sq;
      CUBLAS_CHECK(cublasZdotc(handle, N,
                                s.d_tmp_rot, 1,
                                s.d_tmp_rot, 1,
                                &norm_sq));
      const double beta_k = std::sqrt(cuCreal(norm_sq));
      s.lme[k2 - 1] = beta_k;

      // evec[k2] = d_tmp_rot / beta_k
      if (beta_k > 1.0e-20) {
        const double inv_beta = 1.0 / beta_k;
        CUDA_CHECK(cudaMemcpy(s.d_evec[k2], s.d_tmp_rot, N * CD, D2D));
        CUBLAS_CHECK(cublasZdscal(handle, N, &inv_beta, s.d_evec[k2], 1));
      }

      // Update d_f.
      CUDA_CHECK(cudaMemcpy(s.d_f, s.d_tmp_rot, N * CD, D2D));
    }

    // f) Convergence test.
    //    Re-diagonalize with Nk x Nk block to get Ritz values.
    for (int k = 0; k < Nm; ++k) {
      eval2[k] = s.lmd[k];
      lme2[k]  = s.lme[k];
    }
    s.Qt = Eigen::MatrixXd::Identity(Nm, Nm);
    diagonalize_Eigen(eval2, lme2, Nk, Nm, s.Qt);

    // Allocate a single test vector on the device.
    CuC* d_B;
    CUDA_CHECK(cudaMalloc(&d_B, N * CD));
    CuC* d_HB;
    CUDA_CHECK(cudaMalloc(&d_HB, N * CD));

    int Nconv = 0;
    // Test eigenpairs j=0..Nk-1 (checking all Nk for convergence).
    bool all_conv = true;
    for (int j = 0; j < Nk; ++j) {
      // Construct Ritz vector B = sum_i Qt(i,j) * evec[i] for i=0..Nk-1.
      basisRotateJ(d_B, s.d_evec, s.Qt, j, 0, Nk, Nm, s.d_tmp_rot);

      // Apply HermOp to get H*B.
      HermOp(d_HB, d_B);

      // vnum = Re(<B, H*B>), vden = ||B||^2.
      CuC vnum_c, vden_c;
      CUBLAS_CHECK(cublasZdotc(handle, N, d_B,  1, d_HB, 1, &vnum_c));
      CUBLAS_CHECK(cublasZdotc(handle, N, d_B,  1, d_B,  1, &vden_c));
      const double vnum = cuCreal(vnum_c);
      const double vden = cuCreal(vden_c);
      const double eval_j = vnum / vden;

      // vv = ||H*B - eval_j * B||^2 / evalMaxApprox^2
      // Compute v = H*B - eval_j * B.
      Taxpy_gen<CuC, double, N><<<NBlocks, NThreadsPerBlock>>>(
        d_HB, -eval_j, d_B, d_HB);
      CUDA_CHECK(cudaDeviceSynchronize());

      CuC vv_c;
      CUBLAS_CHECK(cublasZdotc(handle, N, d_HB, 1, d_HB, 1, &vv_c));
      const double vv = cuCreal(vv_c) / (evalMaxApprox * evalMaxApprox);

      bool conv_j = (vv < eresid * eresid);
      if (!conv_j) all_conv = false;
      else ++Nconv;

      std::clog << j << " " << eval_j << " " << 0.0 << " " << eval_j << std::endl;
    }

    CUDA_CHECK(cudaFree(d_B));
    CUDA_CHECK(cudaFree(d_HB));

    std::clog << "# modes converged: " << Nconv << " / " << Nk << std::endl;

    if (all_conv) {
      converged_flag = true;
      break;
    }
  } // end restart loop

  if (!converged_flag) {
    std::cerr << "# IRL did NOT converge after " << MaxIter << " restarts." << std::endl;
  } else {
    std::clog << "# IRL converged at restart " << iter << std::endl;
  }

  // ----------------------------------------------------------
  // Step 4: Final rotation and output.
  // ----------------------------------------------------------
  basisRotate(s.d_evec, s.Qt, 0, Nk, Nm, s.d_tmp_rot, handle);

  // Compute Rayleigh quotients and print.
  CuC* d_B2;
  CUDA_CHECK(cudaMalloc(&d_B2, N * CD));
  CuC* d_HB2;
  CUDA_CHECK(cudaMalloc(&d_HB2, N * CD));

  for (int j = 0; j < Nk; ++j) {
    HermOp(d_HB2, s.d_evec[j]);

    CuC vnum_c, vden_c;
    CUBLAS_CHECK(cublasZdotc(handle, N, s.d_evec[j], 1, d_HB2, 1, &vnum_c));
    CUBLAS_CHECK(cublasZdotc(handle, N, s.d_evec[j], 1, s.d_evec[j], 1, &vden_c));
    const double eval_j = cuCreal(vnum_c) / cuCreal(vden_c);

    s.eval.push_back(eval_j);
    std::clog << j << " " << eval_j << " " << 0.0 << " " << eval_j << std::endl;
  }

  CUDA_CHECK(cudaFree(d_B2));
  CUDA_CHECK(cudaFree(d_HB2));
}
