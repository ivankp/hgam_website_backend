#include "linalg.h"
#include <math.h>

void cholesky(double* A, unsigned N) {
  unsigned col = 0, row = 0;
  for (unsigned k=0; k<N; ++k) {
    double s = 0;
    if (col==row) {
      ++row;
      for (unsigned i=1; i<row; ++i) {
        const double a = A[k-i];
        s += a*a;
      }
      A[k] = sqrt(A[k]-s);
      col = 0;
    } else {
      const unsigned r = col*(col+1u) >> 1u;
      for (unsigned i=0; i<col; ++i)
        s += A[r+i]*A[k+i-col];
      A[k] = (A[k]-s)/A[r+col];
      ++col;
    }
  }
}

void solve_triang(const double* L, double* v, unsigned n) {
  v[0] /= L[0];
  for (unsigned i=1, k=0; i<n; ++i) {
    for (unsigned j=0; j<i; ++j)
      v[i] -= L[++k] * v[j];
    v[i] /= L[++k];
  }
}

void solve_triang_T(const double* L, double* v, unsigned n) {
  unsigned k = n*(n+1u) >> 1u, i = n-1;
  v[i] /= L[--k];
  for (; i; ) {
    --i;
    unsigned k2 = --k;
    for (unsigned j=n-1; j>i; --j) {
      v[i] -= L[k2] * v[j];
      k2 -= j;
    }
    v[i] /= L[k2];
  }
}

void inv_triang(double* L, unsigned n) {
  unsigned i = 0;
  for (unsigned r=0; r<n; ++r) {
    unsigned c2 = 0;
    for (unsigned c=0; c<r; ++c) {
      unsigned j = 0;
      for (; j<c; ++j, ++c2) L[i+j] -= L[i+c]*L[c2];
      L[i+j] *= -L[c2];
      ++c2;
    }
    i += r;
    double l = 1./L[i];
    for (; c2<i; ++c2) L[c2] *= l;
    L[i] = l;
    ++i;
  }
}

double dot(const double* a, const double* b, unsigned n) {
  double x = 0;
  for (unsigned i=0; i<n; ++i)
    x += a[i] * b[i];
  return x;
}

void change_poly_coords(double* c, unsigned n, double a, double b) {
  double p, C;
  for (unsigned i=1; i<n; ++i) {
    C = 1;
    p = b;
    for (unsigned j=1; ; ) {
      C = C * (i+1-j) / j;
      c[i-j] += C * c[i] * p;
      if (j==i) break;
      ++j;
      p *= b;
    }
  }
  p = a;
  for (unsigned i=1; ; ) {
    c[i] *= p;
    if ((++i)>=n) break;
    p *= a;
  }
}

void LT_L(double* L, unsigned n) {
  for (unsigned i=0, k=0; i<n; ++i) {
    for (unsigned j=0; j<=i; ++j, ++k) {
      const unsigned r = i-j;
      unsigned a = k;
      double l = L[a]*L[a+r];
      for (unsigned d=i+1; d<n; ++d) {
        a += d;
        l += L[a]*L[a+r];
      }
      L[k] = l;
    }
  }
}
