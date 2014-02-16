
/** @file hypercluster.h
 * Prototypes for hyperclustering a graph
 * A hyperclustering is a set of overlapping clusters.
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
*/


/** History
 *  2008-09-01: Initial coding 
 */

#include "sparfun.h"
#include <vector>
#include <utility>

int hypercluster_pagerank(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t overlap, size_t minsize, double maxcond,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<double>& conductance,
                          std::vector<int>& order, FILE* statsfile);
                          
int hypercluster_pagerank_full(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t overlap, size_t minsize, double maxcond,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<double>& conductance,
                          std::vector<int>& order, FILE* statsfile);                          
                          
int hypercluster_pagerank_single(sparserow* G, int start, double alpha, 
                                 int max_vol, 
                                 std::vector<int>& cluster);
                                 
int multicluster_pagerank(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t minsize,
                          double maxcond, size_t maxvisits,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<double>& conductance,
                          std::vector<int>& vertices,
                          FILE* statsfile);
