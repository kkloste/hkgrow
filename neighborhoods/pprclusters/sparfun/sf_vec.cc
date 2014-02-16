/** @file sf_vec.cc
 * Routines for working with dense vectors as simple arrays
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
 */

#include "vecfun.h"
#include <assert.h>
#include <math.h>

inline float fabs_t(float x) { return fabsf(x); }
inline double fabs_t(double x) { return fabs(x); }

template <typename ftype, typename stype>
ftype csum_t(const ftype *x, stype n) {
  assert(n >= 0);
  ftype s,e,t,y; s=0.; e=0.;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s,e) shared(x) private(t,y)
  for (stype i = 0; i < n; i++) {
    CSUM2A(x[i],s,e);
  }
  return FCSUM2(s,e);
#else
  for (stype i = 0; i < n; i++) { 
    CSUM2A(x[i],s,e);
  }
  return FCSUM2(s,e);
#endif
}

template <typename ftype, typename stype>
ftype shift_norm_and_norm_diff_t(ftype *y, const ftype *x, const ftype s, 
                               stype n, ftype *pnd, ftype *pny) {
  assert(n >= 0);
  ftype ns=0., ne=0., ds=0., de=0., t, z;
#pragma omp parallel for reduction(+:ns,ne,ds,de) shared(x,y) private(t,z)
  for (stype i = 0; i < n; i++) {
    y[i] += s; CSUM2(fabs_t(y[i]-x[i]),ds,de,t,z); CSUM2(y[i],ns,ne,t,z);
  } 
  if (pnd) { *pnd = FCSUM2(ds,de); }
  if (pny) { *pny = FCSUM2(ns,ne); }
  return FCSUM2(ds,de);
}

float csum(const float *x, int n) { return csum_t(x,n); }
double csum(const double *x, int n) { return csum_t(x,n); }

float shift_norm_and_norm_diff(float *y, const float *x, const float s, 
                               int n, float *pnd, float *pny)
{ return shift_norm_and_norm_diff_t(y,x,s,n,pnd,pny); }  

double shift_norm_and_norm_diff(double *y, const double *x, const double s, 
                               int n, double *pnd, double *pny)
{ return shift_norm_and_norm_diff_t(y,x,s,n,pnd,pny); }

