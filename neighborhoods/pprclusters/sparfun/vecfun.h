/** 
 * @file vecfun.h
 * Prototypes for the vector functions.
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
 */

/** History
 *  2008-07-14: Initial coding
 */

#define ICSUM(s) {s[0]=0.; s[1]=0.;}
#define CSUM(x,s,t,y) { t=s[0]; y=(x)+s[1]; s[0]=t+y; s[1]=(t-s[0])+y; }
#define CSUM2(x,s,e,t,y) { t=s; y=(x)+e; s=t+y; e=(t-s)+y; }
#define FCSUM(s) (s[0]+s[1]) 
#define FCSUM2(s,e) (s+e)
#define CSUMA(x,s) CSUM(x,s,t,y)
#define CSUM2A(x,s,e) CSUM2(x,s,e,t,y)

double csum(const double *x, int n);
float csum(const float *x, int n);

double shift_norm_and_norm_diff(double *y, const double *x, const double s, 
                               int n, double *pnd, double *pny);

