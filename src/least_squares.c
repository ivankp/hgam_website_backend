#include "least_squares.h"
#include "linalg.h"
#include <string.h>

void linear_least_squares(
  unsigned nx, // number of measured values
  unsigned np, // number of parameters
  const double* A, // design matrix, [np][nx]
  const double* y, // observed values, [nx]
  const double* u, // variances, [nx] or null
  double* workspace, // N = np*(np+1)/2, [N ( + nx if u )]
  double* P, // fitted coefficients (parameters), [np]
  double* cov // covariance matrix, [N]
) {
  const unsigned N = np*(np+1u) >> 1u; // upper triangular number
  double* const L = workspace;
  double* const V = L + N;

  if (!u) { // no input variances

    // At A
    for (unsigned i=0, r1=0, k=0; i<np; ++i, r1+=nx) {
      for (unsigned j=0, r2=0; j<=i; ++j, r2+=nx, ++k) {
        double l = 0;
        for (unsigned x=0; x<nx; ++x)
          l += A[r1+x] * A[r2+x];
        L[k] = l;
      }
    }

    // (At A)^-1
    cholesky(L,N);

    // At y
    for (unsigned i=0, row=0; i<np; ++i, row+=nx) {
      double p = 0;
      for (unsigned x=0; x<nx; ++x)
        p += A[row+x] * y[x];
      P[i] = p;
    }

  } else { // input variances provided

    // V^-1
    for (unsigned x=0; x<nx; ++x)
      V[x] = 1./u[x];

    // At V^-1 A
    for (unsigned i=0, r1=0, k=0; i<np; ++i, r1+=nx) {
      for (unsigned j=0, r2=0; j<=i; ++j, r2+=nx, ++k) {
        double l = 0;
        for (unsigned x=0; x<nx; ++x)
          l += A[r1+x] * A[r2+x] * V[x];
        L[k] = l;
      }
    }

    // (At V^-1 A)^-1
    cholesky(L,N);

    // V^-1 y
    for (unsigned x=0; x<nx; ++x)
      V[x] *= y[x];

    // At V^-1 y
    for (unsigned i=0, row=0; i<np; ++i, row+=nx) {
      double p = 0;
      for (unsigned x=0; x<nx; ++x)
        p += A[row+x] * V[x];
      P[i] = p;
    }

  }

  solve_triang  (L,P,np); // solve p = L^-1 p
  solve_triang_T(L,P,np); // solve p = LT^-1 p

  // covariance matrix
  // 0 1 3 6   00 10 20 30
  //   2 4 7      11 21 31
  //     5 8         22 32
  //       9            33

  if (cov) {
    inv_triang(L,np); // invert
    LT_L(L,np); // multiply LT by L
    memcpy(cov,L,N*sizeof(double));
  }
}

double linear_least_squares_chi2(
  unsigned nx, // number of measured values
  unsigned np, // number of parameters
  const double* A, // design matrix, [np][nx]
  const double* y, // observed values, [nx]
  const double* u, // variances, [nx] or null
  const double* P  // fitted coefficients (parameters), [np]
) {
  double chi2 = 0;
  for (unsigned i=0; i<nx; ++i) {
    double w = y[i];
    for (unsigned j=0; j<np; ++j) {
      w -= A[nx*j+i] * P[j];
    }
    w *= w;
    if (u) w /= u[i];
    chi2 += w;
  }
  return chi2;
}
