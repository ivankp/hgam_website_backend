#ifndef IVAN_LINALG_H
#define IVAN_LINALG_H

void cholesky(double* A, unsigned N);
void solve_triang(const double* L, double* v, unsigned n);
void solve_triang_T(const double* L, double* v, unsigned n);
void inv_triang(double* L, unsigned n);
double dot(const double* a, const double* b, unsigned n);

void LT_L(double* L, unsigned n);

void change_poly_coords(double* c, unsigned n, double a, double b);

#endif
