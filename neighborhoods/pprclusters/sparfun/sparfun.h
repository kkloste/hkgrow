/** 
 * @file sparfun.h
 * Prototypes for the sparfun library functions.
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
 */

#ifndef SPARFUN_SPARFUN_H
#define SPARFUN_SPARFUN_H

/** History
 *  2008-07-03: Initial coding
 */
#include <stdio.h>

struct sparserow {
  int n, m; // number of 
  int *ai, *aj;
  double *a;
};

int load_bsmat(const char* filename, sparserow* s, bool load_vals=true);
int load_smat(const char* filename, sparserow* s, bool load_vals=true);

int load_generic(const char* filename, sparserow* s, bool load_vals=true);

/** Return the degree of a node 
 */
inline int sr_degree(const sparserow *s, int u) {
  return s->ai[u+1]-s->ai[u];
}
inline int sr_nedges(const sparserow *s) {
  return s->ai[s->m];
}
int sr_graph_volume(const sparserow *s, int *verts, size_t n);
int sr_graph_cutsize(const sparserow *s, int *verts, size_t n);
int sr_graph_core_numbers(const sparserow *s, int *core_map);

sparserow* sf_sym_or(sparserow *s);
int sf_pack(sparserow *s);
int sf_sortvals(sparserow *r);
sparserow* sf_transpose(sparserow *s, bool values=true);
sparserow* sf_row_subset(sparserow *r, size_t *rows, size_t nrows);
sparserow* sf_index_int_vector_elements(int* vector, size_t nelem, 
                size_t maxvalue);
                
sparserow* sf_rowcol_subset(sparserow *r, size_t *set, size_t nset);                

void free_sparse(sparserow* s);
int sparserow_alloc(sparserow *s, int n, int nz, bool vals);
sparserow* sparserow_fullalloc(int n, int nz, bool vals);

#endif /* SPARFUN_SPARFUN_H */
