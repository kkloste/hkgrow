/** @file sf_matrix.cc
 * Routines for working with sparserow types as matrices
 */

/*
 * David F. Gleich
 * Copyright, University of British Columbia, 2010
 * 
 * History
 * -------
 * :2010-01-27: Added sf_row_subset
 */

#include "sparfun.h"
#include <assert.h>
#include <stdlib.h> 
#include <string.h> // used for memset
#include <algorithm> // used for std::sort
#include <vector>

/** Create a sparse matrix to efficient access elements of a vector by value.
 * This only applies when the values of a vector are integers and there
 * are only a few of them.  We just build a sparse matrix to do this...
 * 
 * Example
 * -------
 * vector = [0 1 2 0 0 1 1 1 1 ]
 * index = [ 0 => 0, 3, 4
 *           1 => 1, 5, 6, 7, 8
 *           2 => 2 ]
 * 
 * @param vector the vector that we are going to index
 * @param nelem the length of the vector 
 * @param maxvalue the maximum element in the vector.
 */
sparserow* sf_index_int_vector_elements(int* vector, size_t nelem, 
                size_t maxvalue)
{
    size_t i,j,nzi;
    sparserow *s = sparserow_fullalloc(maxvalue, nelem, false);
    if (s==NULL) { return NULL; } // TODO add error
    memset(s->ai, 0, sizeof(int)*(s->m+1));
    for (i=0; i<nelem; i++) {
        // we have a non-zero vector[i], i
        assert(vector[i]>=0);
        assert(vector[i]<maxvalue);
        s->ai[vector[i]+1]++;
    }
    for (nzi=0, j=0; j<s->m+1; j++) { s->ai[j] = (nzi += s->ai[j]); }
    for (i=0; i<nelem; i++) {
        j = s->ai[vector[i]];
        s->aj[s->ai[vector[i]]] = (int)i;
        s->ai[vector[i]]++;
    }
    for (j=s->m; j>0; j--) { s->ai[j] = s->ai[j-1]; }
    s->ai[0] = 0;
    return s;
}
        

sparserow* sf_sym_or(sparserow *r) 
{
  int i,j,nzi;
  sparserow *s=sparserow_fullalloc(r->m>r->n?r->m:r->n ,2*r->ai[r->m],false);
  if (s == NULL) { return NULL; } // TODO add error
  s->n = s->m; // assign the largest size
  memset(s->ai, 0, sizeof(int)*(s->m+1));
  for (i=0; i<r->m; ++i) { 
    for (nzi=r->ai[i];nzi<r->ai[i+1];++nzi) { s->ai[i+1]++; s->ai[r->aj[nzi]+1]++; }
  }
  for (nzi=0, j=0; j<s->m+1; j++) { s->ai[j] = (nzi += s->ai[j]); }
  for (i=0; i<r->m; ++i) { 
    for (nzi=r->ai[i];nzi<r->ai[i+1];++nzi) { 
      j=r->aj[nzi]; 
      s->aj[s->ai[i]] = j; s->ai[i]++;
      s->aj[s->ai[j]] = i; s->ai[j]++;
    }
  }
  for (j=s->m; j>0; j--) { s->ai[j] = s->ai[j-1]; }
  s->ai[0] = 0;
  return s;
}

sparserow* sf_sym_max(sparserow *r)
{
  if (r->a == NULL) {
    return sf_sym_or(r);
  } else {
    return NULL;
  }
}

/** Compute a new sparse row matrix from a subset of rows.
 * @param r the initial matrix
 * @param rows the subset of rows, where rows[i] is the index of a row
 * @param nrows the number of total rows in the subset
 * @return a new sparse matrix only over these rows
 */
sparserow* sf_row_subset(sparserow *r, size_t *rows, size_t nrows)
{
    // count non-zeros in the new matrix
    size_t newnz = 0;
    bool vals = false;
    if (r->a) { vals = true; }
    for (size_t i=0; i<nrows; i++) {
        size_t ri = rows[i];
        newnz += r->ai[ri+1] - r->ai[ri];
    }
    sparserow* s = sparserow_fullalloc(nrows, newnz, vals);
    if (s==NULL) { return NULL; } // TODO Add error
    s->n = r->n;
    size_t curnz = 0;
    // copy non-zeros
    for (size_t i=0; i<nrows; i++) {
        size_t ri = rows[i];
        s->ai[i] = curnz;
        for (int nzi = r->ai[ri]; nzi < r->ai[ri+1]; ++nzi) {
            s->aj[curnz] = r->aj[nzi];
            if (vals) { s->a[curnz] = r->a[nzi]; }
            curnz++;
        }
    }
    s->ai[nrows] = curnz;
    return s;
}

/** Compute a new sparse matrix from a subset of rows and columns
 * 
 * The order of the new matrix is the same as the order of the rows
 * in set.
 * 
 * @param r the initial matrix
 * @param set the subset of rows and columns 
 * @param nset the size of the set
 */
sparserow* sf_rowcol_subset(sparserow *r, size_t *set, size_t nset)
{
  bool values = r->a != NULL;
  // setup iwork array to hold new indices
  int* iwork = (int*)malloc(sizeof(int)*r->n);
  if (!iwork) { return NULL; }  // TODO add error
  
  // build map for fast lookup 
  memset(iwork, -1, sizeof(int)*r->n);
  for (size_t i=0; i<nset; ++i) {
    assert(set[i] < r->n);
    assert(set[i] < r->m);
    iwork[set[i]] = i;
  }
  
  // compute the non-zeros of the new matrix
  size_t newnz = 0;
  for (size_t i=0; i<nset; ++i) {
    size_t ri = set[i];
    for (int nzi=r->ai[ri]; nzi<r->ai[ri+1]; ++nzi) {
      if (iwork[r->aj[nzi]] >= 0) {
        newnz += 1;
      }
    }
  }
  
  // allocate memory
  sparserow *s=sparserow_fullalloc(nset,newnz,values);
  if (s == NULL) { free(iwork); return NULL; } // TODO add error
  s->n = nset;
  
  size_t curnz = 0;
  // copy non-zeros
  for (size_t i=0; i<nset; i++) {
    size_t ri = set[i];
    s->ai[i] = curnz;
    for (int nzi = r->ai[ri]; nzi < r->ai[ri+1]; ++nzi) {
      if (iwork[r->aj[nzi]] >= 0) {
        s->aj[curnz] = iwork[r->aj[nzi]];
        if (values) { s->a[curnz] = r->a[nzi]; }
        curnz++;
      }
    }
  }
  s->ai[nset] = curnz;
  
  assert(curnz == newnz);
  
  free(iwork);
  return s;
}

template <typename T>
class sortvals_functor {
public:
  T* _off;
  void offset(T *off) { _off = off; }
  bool operator() (const int a, const int b) const { return _off[a]<_off[b];  }
};

/** Sort the elements of the row by their value
 * todo: add an option to sort increasing or decreasing
 */
int sf_sortvals(sparserow *r) 
{
  if (!r->a) { return 1; } // TODO return error?
  std::vector<int> perm;
  std::vector<int> jcopy;
  std::vector<double> acopy;
  sortvals_functor<double> f;
  for (int i=0; i<r->m; i++) {
    // ugh, this is so ridiculous that nothing better is easier
    // std::sort should work on the zipped iterator, but the language
    // semantics don't allow it.  So we resort to this uglyness.
    int nedges = r->ai[i+1]-r->ai[i];
    if (perm.size()<nedges) { 
      perm.resize(nedges); jcopy.resize(nedges); acopy.resize(nedges); 
    }
    std::copy(&r->aj[r->ai[i]],&r->aj[r->ai[i]]+nedges,jcopy.begin());
    std::copy(&r->a[r->ai[i]],&r->a[r->ai[i]]+nedges,acopy.begin());
    f.offset(&acopy[0]);
    for (int ei=0; ei<nedges; ++ei) {
      perm[ei]=ei;
    }
    std::sort(&perm[0], &perm[0]+nedges, f);
    // apply perm
    for (int ei=0; ei<nedges; ++ei) {
      r->aj[r->ai[i]+ei] = jcopy[perm[ei]];
      r->a[r->ai[i]+ei] = acopy[perm[ei]];
    }
  }
  return (0);
}

/** Sort the elements of the row by their value
 * todo: add an option to sort increasing or decreasing
 */
int sf_sortrows(sparserow *r) 
{
  if (r->a) { return 1; } // TODO return error?
  std::vector<int> perm;
  std::vector<int> jcopy;
  std::vector<double> acopy;
  sortvals_functor<int> f;
  for (int i=0; i<r->m; i++) {
    // ugh, this is so ridiculous that nothing better is easier
    // std::sort should work on the zipped iterator, but the language
    // semantics don't allow it.  So we resort to this uglyness.
    int nedges = r->ai[i+1]-r->ai[i];
    if (perm.size()<nedges) { 
      perm.resize(nedges); jcopy.resize(nedges); acopy.resize(nedges); 
    }
    std::copy(&r->aj[r->ai[i]],&r->aj[r->ai[i]]+nedges,jcopy.begin());
    std::copy(&r->a[r->ai[i]],&r->a[r->ai[i]]+nedges,acopy.begin());
    f.offset(&jcopy[0]);
    for (int ei=0; ei<nedges; ++ei) {
      perm[ei]=ei;
    }
    std::sort(&perm[0], &perm[0]+nedges, f);
    // apply perm
    for (int ei=0; ei<nedges; ++ei) {
      r->aj[r->ai[i]+ei] = jcopy[perm[ei]];
      r->a[r->ai[i]+ei] = acopy[perm[ei]];
    }
  }
  return (0);
}

/** Remove duplicate non-zeros 
 */
int sf_pack(sparserow* r) 
{
  assert(r->a==NULL); // TODO remove this assert when we handle values
  int* iwork = (int*)malloc(sizeof(int)*r->n);
  if (!iwork) { return -1; }  // TODO add error
  memset(iwork,-1,sizeof(int)*r->n);
  int i,j,nzi,curnz;
  curnz=0;
  for (i=0;i<r->m;++i) {
    int rstart=curnz;
    for (nzi=r->ai[i];nzi<r->ai[i+1];++nzi) {
      j=r->aj[nzi];
      if (iwork[j]==-1) { // new column
        iwork[j]=curnz; r->aj[curnz]=j; curnz++;
      }
    }
    r->ai[i]=rstart;
    for (nzi=r->ai[i];nzi<curnz;++nzi) { iwork[r->aj[nzi]]=-1; }
  }
  r->ai[r->m]=curnz;
  free(iwork);
  return 0;
}

/** Tranpose a sparse matrix
 */
sparserow* sf_transpose(sparserow *r, bool values) 
{
  int i,j,nzi;
  values = values && (r->a != NULL);
  
  // allocate memory
  sparserow *s=sparserow_fullalloc(r->n,r->ai[r->m],values);
  if (s == NULL) { return NULL; } // TODO add error
  
  // set ai to 0, then compute row counts for s from column counts for r
  memset(s->ai, 0, sizeof(int)*(s->m+1));
  for (i=0; i<r->m; ++i) { 
    for (nzi=r->ai[i];nzi<r->ai[i+1];++nzi) { 
      s->ai[r->aj[nzi]+1]++; // the +1 is to make the cum-sum easier
    }
  }
  // cum-sum the vector s->ai, this produces the starting location for all nz
  for (nzi=0, j=0; j<s->m+1; j++) { s->ai[j] = (nzi += s->ai[j]); }
  for (i=0; i<r->m; ++i) { 
    for (nzi=r->ai[i];nzi<r->ai[i+1];++nzi) { 
      j=r->aj[nzi]; 
      s->aj[s->ai[j]] = i; 
      if (values) { s->a[s->ai[j]] = r->a[nzi]; }
      s->ai[j]++;
    }
  }
  // return the vectorto its prior state
  for (j=s->m; j>0; j--) { s->ai[j] = s->ai[j-1]; }
  s->ai[0] = 0;
  return s;
}
