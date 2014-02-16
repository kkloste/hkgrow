#!/usr/bin/env python

"""
dgg.py
======

dgg or dgleich-graph is a quick wrapper to interface igraph with my standard
graph formats.
"""

import igraph
import time

class Timer:
   def __enter__(self): self.start = time.time()
   def __exit__(self, *args): print time.time() - self.start


def read_smat(filename,set_undirected=True):
    """ Load a graph in SMAT format.
    
    G = dgg.read_smat(filename,set_undirected=True)
    loads filename as an igraph graph based on the smat format.
    
    
    @param filename the path to the file
    @param set_undirected a flag to check if the graph is undirected and set igraph
      accordingly.  If this is false, the output is always a directed igraph graph.
    
    @return the igraph graph.  
    """
    
    f = open(filename,'rU') # use universal newlines
    hdr = f.readline();
    parts = hdr.split();
    nverts = int(parts[0])
    ncols = int(parts[1])
    nnz = int(parts[2])
    if nverts != ncols:
        raise ValueError(
            'read_smat needs nrows (%i) = ncols (%i) in smat for graph'%(
            nverts, ncols))
    edgelist = []
    for line in f:
        parts = line.split()
        edgelist.append((int(parts[0]),int(parts[1])))
        
    G = igraph.Graph(nverts,edges=edgelist,directed=True)
    if G.ecount() != nnz:
        raise ValueError(
            'graph has %i edges but smat file has %i nnz'%(
            G.ecount(), nnz))
            
    if G.vcount() > nverts:
        raise ValueError(
            'graph has more vertices (%i) than smat file has rows %i'%(
            G.vcount(), nverts))        
    
    if set_undirected:
        dir = False
        for medge in G.is_mutual():
            if medge is False:
                dir = True
                break
        if not dir:
            G.to_undirected()
    
    if G.vcount() < nverts:
        G.add_vertices(nverts-G.vcout())

    return G
    
def cond(G,set):
    Gvol = sum(G.degree(type=igraph.OUT))
    Svol = vol(G,set)
    return float(cutsize(G,set))/float(min(Svol,Gvol-Svol))
    
def cut(G,S):
    fastset = set(S)
    cut = []
    for u in fastset:
        for v in G.neighbors(u,type=igraph.OUT):
            if v not in fastset:
                cut.append((u,v))
    return cut
    
def cutsize(G,S):
    fastset = set(S)
    cut = 0
    for u in fastset:
        for v in G.neighbors(u,type=igraph.OUT):
            if v not in fastset:
                cut += 1
    return cut
    
def vol(G,set):
    degs = G.degrees(set,type=igraph.OUT)
    return sum(degs)
    
    
def clustercoeffs(G):
    ccs = [0. for v in xrange(G.vcount())]
    marker = [False for v in xrange(G.vcount())]
    for v in xrange(G.vcount()):
        ntris = 0
        
        d = 0
        for u in G.neighbors(v,type=igraph.OUT):
            marker[u] = True
            if u != v:
                d += 1
            
        for u in G.neighbors(v,type=igraph.OUT):
            for w in G.neighbors(u,type=igraph.OUT):
                if marker[w] and u != w and w != v:
                    ntris += 1
                    
        for u in G.neighbors(v,type=igraph.OUT):
            marker[u] = False
        
        if d>1:
            ccs[v] = (float(ntris)/(float(d)*float(d-1)))
        else:
            ccs[v] = 0.
    return ccs
        