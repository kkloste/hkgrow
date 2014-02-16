/** @file sf_graph.cc
 * Routines for working with sparserow types as graphs
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
 */

#ifdef __APPLE__
#include <tr1/unordered_set>
#else
#include <unordered_set>
#endif 

#if defined(_WIN32) || defined(_WIN64) || defined(__APPLE__)
#define tr1ns std::tr1
#else
#define tr1ns std
#endif

#include <vector>
#include "sparfun.h"
#include "sparvec.h"  // include macros to make working with tr1 better

/** Compute the volume of a set of vertices.
 * 
 * For a given set of vertices, the volume is the sum of their
 * degrees.  This function will compute that quantity for
 * an array of vertices.
 * 
 * @param c the graph
 * @param verts an array of vertices numbers in the set
 * @param nverts the number of vertices in the set (size of verts)
 * @return the volume of the set
 */ 
int sr_graph_volume(const sparserow* c, int* verts, size_t nverts) {
    int vol = 0;
    for (size_t i=0; i<nverts; i++) {
        int v = verts[i];
        vol += c->ai[v+1]-c->ai[v];
    }
    return vol;
}

/** Compute the cut-size (boundary-size) of a set of vertices.
 * 
 * For a given set of vertices, the volumen is the sum of their
 * degrees.  This function will compute that quantity for
 * an array of vertices.
 * 
 * @param c the graph
 * @param verts an array of vertices numbers in the set
 * @param nverts the number of vertices in the set (size of verts)
 */ 
int sr_graph_cutsize(const sparserow* c, int* verts, size_t nverts) {
    int cut = 0;
    tr1ns::unordered_set<int> verts_set;
    for (size_t i=0; i<nverts; i++) {
        verts_set.insert(verts[i]);
    }
    for (size_t i=0; i<nverts; i++) {
        // look at all neighbors
        int v = verts[i];
        for (int nzi=c->ai[v]; nzi<c->ai[v+1]; ++nzi) {
            int neighbor = c->aj[nzi];
            if (verts_set.count(neighbor) == 0) {
                cut += 1;
            }
        }
    }
    return cut;
}

/**
 * @param g the sparserow graph
 * @param core_map an array of size g->n to hold the core numbers
 * @return the largest core number
 */
int sr_graph_core_numbers(const sparserow* g, int *core_map) 
{
    // compute the maximum degree and store degrees is core_map
    int max_degree = 0;
    for (int i=0; i<g->n; ++i) {
        core_map[i] = sr_degree(g,i);
        if (core_map[i] > max_degree) {
            max_degree = core_map[i];
        }
    }
    
    // in what follows, we are bucket sorting vertices by degree
    
    // store the vertices in bins by their degree
    // allocate two extra locations to ease boundary cases
    std::vector<int> bin(max_degree+2);
    for (int i=0; i<g->n; ++i) {
        bin[core_map[i]]++; // count the number of vertices with degree i
    }
    // compute the cum-sum of the degree bins to do the bucket sort placement
    int cur_pos = 0;
    for (int cur=0; cur<max_degree+1; ++cur) {
        int tmp = bin[cur];
        bin[cur] = cur_pos;
        cur_pos += tmp;
    }
    // complete the bucket sort by storing vertex ids in sorted order
    std::vector<int> vert(g->n);
    std::vector<int> pos(g->n);
    for (int i=0; i<g->n; ++i) {
        int p = bin[core_map[i]];
        pos[i] = p;
        vert[p] = i;
        ++bin[core_map[i]];
    }

    // now fix bin to restore its original contents
    for (int cur=max_degree+1; cur>0; --cur) {
        bin[cur] = bin[cur-1];
    }
    bin[0] = 0;

    // now run through and simulate removing the vertices
    for (int i=0; i<g->n; ++i) {
        // remove vertex vert[i];
        int v = vert[i];
        int v_cn = core_map[v];
        for (int nzi=g->ai[v]; nzi<g->ai[v+1]; ++nzi) {
            int u = g->aj[nzi];
            if (core_map[u] > v_cn) {
                // u is still in the graph
                int u_cn = core_map[u];
                int u_pos = pos[u];
                
                // now find the first vertex with the same core_number as u
                // (this resorts)
                int w_pos = bin[u_cn];
                int w = vert[w_pos];
                if (u!=v) {
                    // swap u and w
                    pos[u]=w_pos;
                    vert[w_pos] = u;
                    pos[w]=u_pos;
                    vert[u_pos] = w;              
                }
                ++bin[u_cn];
                --core_map[u];
            }
        }
    }
    return core_map[g->n-1];
}
