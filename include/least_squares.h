#ifndef IVAN_LEAST_SQUARES_H
#define IVAN_LEAST_SQUARES_H

#ifdef __cplusplus
extern "C" {
#endif

void linear_least_squares(
  unsigned nx, // number of measured values
  unsigned np, // number of parameters
  const double* A, // design matrix, [np][nx]
  const double* y, // observed values, [nx]
  const double* u, // variances, [nx] or null
  double* workspace, // N = np*(np+1)/2, [N ( + nx if u )]
  double* P, // fitted coefficients (parameters), [np]
  double* cov // covariance matrix, [N]
);

double linear_least_squares_chi2(
  unsigned nx, // number of measured values
  unsigned np, // number of parameters
  const double* A, // design matrix, [np][nx]
  const double* y, // observed values, [nx]
  const double* u, // variances, [nx] or null
  const double* P  // fitted coefficients (parameters), [np]
);

#ifdef __cplusplus
}
#endif

#endif
