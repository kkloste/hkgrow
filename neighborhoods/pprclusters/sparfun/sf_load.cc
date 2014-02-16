/**
 * @file sf_load.cc
 * Sparse matrix load routines
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
 */

/** History
 *  2008-07-03: Initial coding
 *  2008-09-01: Changed to use fastio routines
 */

#include "sparfun.h"
#include "fastio.h"

#include <stdlib.h> 
#pragma warning(disable:4996)
#include <stdio.h>

#include <assert.h>
#include <string>

void free_sparse(sparserow *s) {
  if (!s) { return; }
  if (s->a) { free(s->a); s->a = NULL; }
  if (s->ai) { free(s->ai); s->ai = NULL; }
  if (s->aj) { free(s->aj); s->aj = NULL; }
  s->m = 0;
  s->n = 0;
}

sparserow* sparserow_fullalloc(int n, int nz, bool vals) {
  sparserow *s = (sparserow*)malloc(sizeof(sparserow));
  if (!s) { return (NULL); }
  assert(n>=0); assert(nz >= 0);

  // cannot be replaced by sparserow_alloc because of the
  // the free in sparserow_alloc
  s->ai = (int*)malloc(sizeof(int)*(n+1));
  s->aj = (int*)malloc(sizeof(int)*(nz));
  if (vals) { s->a = (double*)malloc(sizeof(double)*nz); }
  else { s->a = NULL; }

  s->m=n;

  return (s);
}


int sparserow_alloc(sparserow *s, int n, int nz, bool vals) {
  free_sparse(s);
  assert(s);
  assert(n>=0); assert(nz >= 0);

  s->ai = (int*)malloc(sizeof(int)*(n+1));
  s->aj = (int*)malloc(sizeof(int)*(nz));
  if (vals) { s->a = (double*)malloc(sizeof(double)*nz); }
  else { s->a = NULL; }

  return (0);
}

int sf_idone(int rval, sparserow *s) 
{
  if (s != NULL) { free_sparse(s); }
  return (rval);
}
int sf_ifulldone(int rval, sparserow* s) {
  if (s != NULL) { free_sparse(s); free(s); }
  return (rval);
}  

/** Load a sparse matrix in bsmat format.
 * 
 * @return 0 on success
 */
int load_bsmat(const char *filename, sparserow* s, bool load_vals)
{
  buffered_reader r(filename, 1<<16);

  int m, n, nz;
  size_t rb; const size_t skip = sizeof(int)+sizeof(double);
  if (r.read((char*)&m, sizeof(int)) &&
      r.read((char*)&n, sizeof(int)) &&
      r.read((char*)&nz, sizeof(int))) {
    if (sparserow_alloc(s, m, nz, load_vals)==0) {
      s->m = m; 
      s->n = n;
      memset(s->ai, 0, sizeof(int)*(s->m+1));
      int i[2], j, nzi = nz;
      while (nzi-- > 0) {
        rb = r.read((char*)&i, sizeof(int), 1); 
        if (rb!=1) return sf_idone(-1,s);
        if (r.seekf(skip) != 0) { return sf_idone(-1,s); }
        if (i[0] >= 0 && i[0] < m) { s->ai[i[0]+1]++; }
        else { return sf_idone(-10,s); }
      }
      for (nzi=0, j=0; j<s->m+1; j++) {
        s->ai[j] = (nzi += s->ai[j]);
      }
      if (r.seek(3*sizeof(int), SEEK_SET) != 0) { return sf_idone(-1,s); }
      nzi = nz; double a;
      while (nzi-- > 0) {
        rb = r.read(&i, sizeof(int), 2); if (rb!=2) return sf_idone(-1,s);
        rb = r.read(&a, sizeof(double), 1); if (rb!=1) return sf_idone(-1,s);
        if (i[0] >= 0 && i[0] < m && i[1] >= 0 && i[1] < n) {
          s->aj[s->ai[i[0]]] = i[1];
        } else { return sf_idone(-10,s); }
        if (load_vals) {s->a[s->ai[i[0]]] = a;}
        s->ai[i[0]]++;
      }

      for (j=s->m; j>0; j--) {
        s->ai[j] = s->ai[j-1];
      }
      s->ai[0] = 0;

      return (0);
    } else {
      return (-3);
    }
  } else {
    return (-1);
  }
}

/** Load a sparse matrix in bsmat format.
 * 
 * @return 0 on success
 */
int load_smat(const char *filename, sparserow* s, bool load_vals)
{
  FILE *f = fopen(filename, "rt");
  if (!f) {
    return (-1);
  }

  int m, n, nz;
  if (fscanf(f, "%u %u %u", &m, &n, &nz)==3) {
    if (sparserow_alloc(s, m, nz, load_vals)==0) {
      s->m = m; 
      s->n = n;
      memset(s->ai, 0, sizeof(int)*(s->m+1));
      int i, j, nzi = nz;
      double a;
      while (nzi-- > 0) {
        if (fscanf(f, "%u %u %lf", &i, &j, &a) != 3) {return sf_idone(-1,s);}
        if (i >= 0 && i < m) { s->ai[i+1]++; }
        else { return sf_idone(-10,s); }
      }
      for (nzi=0, j=0; j<s->m+1; j++) {
        s->ai[j] = (nzi += s->ai[j]);
      }
      if (fseek(f, 0, SEEK_SET) != 0) { return sf_idone(-1,s); }
      if (fscanf(f, "%u %u %u", &m, &n, &nz)!=3) { return sf_idone(-1,s); }
      nzi = nz; 
      while (nzi-- > 0) {
        if (fscanf(f, "%u %u %lf", &i, &j, &a) != 3) {return sf_idone(-1,s);}
        if (i >= 0 && i < m && j >= 0 && j < n) {
          s->aj[s->ai[i]] = j;
        } else { return sf_idone(-10,s); }
        if (load_vals) {s->a[s->ai[i]] = a;}
        s->ai[i]++;
      }
      for (j=s->m; j>0; j--) {
        s->ai[j] = s->ai[j-1];
      }
      s->ai[0] = 0;

      return (0);
    } else {
      return (-3);
    }
  } else {
    return (-1);
  }
}


/** Load a file based on its filetype
 */
int load_generic(const char* filename, sparserow* s, bool load_vals)
{
  std::string sfilename(filename);
  std::string::size_type lastdot = sfilename.rfind('.');
  if (lastdot != std::string::npos) {
    std::string ext = sfilename.substr(lastdot+1);
    if (ext.compare("smat")==0) {
      return load_smat(filename, s, load_vals);
    } else if (ext.compare("smat-sym")==0) {
      return load_smat(filename, s, load_vals);
    } else if (ext.compare("smat-cc")==0) {
      return load_smat(filename, s, load_vals);
    } else if (ext.compare("smat-sym-cc")==0) {
      return load_smat(filename, s, load_vals);
    } else if (ext.compare("eg2")==0) {
      return load_smat(filename, s, load_vals);
    } else if (ext.compare("bsmat")==0) {
      return load_bsmat(filename, s, load_vals);
    }
  }
  // getting here indicates an error
  return (1);
}

