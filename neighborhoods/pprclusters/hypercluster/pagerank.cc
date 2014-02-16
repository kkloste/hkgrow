/** @file pagerank.cc
 * Prototypes for hyperclustering a graph using PageRank partitions
 */

/*
 * David F. Gleich
 * Copyright, Microsoft Corporation, 2008
*/


/** 
 * History
 * -------
 * :2008-09-01: Initial coding 
 * :2010-02-02: Fixed bug with target_vol overflowing an int.  (Happened
 * with 100 sequential trials < 0.1!
 * :2010-07-28: Added arbitrary examination order
 * :2011-10-13: Fixed array_queue
 */

#define NOMINMAX

#include <iostream>

#include "hypercluster.h"
#include "sparvec.h"
#include "sparfun_util.h"
#include <queue>
#include <assert.h>
#include <limits>
#include <algorithm>
#include <stdio.h>
#include <math.h>


/** A replacement for std::queue<int> using a circular buffer array */
class array_queue {
    public:
    std::vector<int> array;
    size_t max_size;
    size_t head, tail;
    size_t cursize;
    array_queue(size_t _max_size)
    : max_size(_max_size), array(_max_size), head(0), tail(0), cursize(0)
    {}
    
    void empty() {
        head = 0;
        tail = 0;
        cursize = 0;
    }
    
    size_t size() {
        return cursize;
    }
    
    void push(int i) {
        assert(size() < max_size);
        array[tail] = i;
        tail ++;
        if (tail == max_size) {
            tail = 0;
        }
        cursize ++;
    }
    
    int front() {
        assert(size() > 0);
        return array[head];
    }
    
    void pop() {
        assert(size() > 0);
        head ++;
        if (head == max_size) {
            head = 0;
        }
        cursize --;
    }
};

template <class Queue>
int compute_local_pagerank(sparserow *s, sparsevec& r, sparsevec& p, 
    double alpha, double epsilon, int max_push_count, Queue& q) 
{
  for (sparsevec::map_type::iterator it=r.map.begin(),itend=r.map.end();
        it!=itend;++it){
    if (it->second > epsilon*sr_degree(s,it->first)) {
      q.push(it->first);
    }
  }

  int push_count = 0;
  while (q.size()>0 && push_count < max_push_count) {
    push_count += 1;
    int u = q.front();
    q.pop();
    int du = sr_degree(s, u);
    double moving_probability = r.map[u] - 0.5*epsilon*(double)du;
    r.map[u] = 0.5*epsilon*(double)du;
    p.map[u] += (1.-alpha)*moving_probability;

    double neighbor_update = alpha*moving_probability/(double)du;

    for (int nzi=s->ai[u]; nzi<s->ai[u+1]; nzi++) {
      int x = s->aj[nzi];
      int dx = sr_degree(s, x);
      double rxold = r.get(x);
      double rxnew = rxold + neighbor_update;
      r.map[x] = rxnew;
      if (rxnew > epsilon*dx && rxold <= epsilon*dx) {
        q.push(x);
      }
    }
  }
  
  return (push_count);
}

/**
 * Todo: document this function (sigh...)
 * Add statistics on what it's doing...
 * @return the number of steps
 */
int compute_local_pagerank(sparserow *s, sparsevec& r, sparsevec& p, 
    double alpha, double epsilon, int max_push_count) 
{
  std::queue<int> q;
  return compute_local_pagerank(s, r, p, alpha, epsilon, max_push_count, q);
}

struct greater2nd {
  template <typename P> bool operator() (const P& p1, const P& p2) {
    return p1.second > p2.second;
  }
};

void cluster_from_sweep(sparserow* G, sparsevec& p, int target_vol, 
      std::vector<int>& cluster, double *outcond, int* outvolume)
{
  // now we have to do the sweep over p in sorted order by value
  typedef std::vector< std::pair<int, double> > vertex_prob_type;
  vertex_prob_type prpairs(p.map.begin(), p.map.end());
  std::sort(prpairs.begin(), prpairs.end(), greater2nd());

  // compute cutsize, volume, and conductance
  std::vector<double> conductance(prpairs.size());
  std::vector<int> volume(prpairs.size());

  size_t i=0;
  tr1ns::unordered_map<int,size_t> rank;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    rank[it->first] = i;
  }
  //printf("support=%i\n",prpairs.size());
  int total_degree = G->ai[G->m];
  int curcutsize = 0;
  int curvolume = 0;
  i=0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    int v = it->first;
    int deg = G->ai[v+1]-G->ai[v];
    int change = deg;
    for (int nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
      int nbr = G->aj[nzi];
      if (rank.count(nbr) > 0) {
        if (rank[nbr] < rank[v]) {
          change -= 2;
        }
      }
    }
    curcutsize += change;
    if (curvolume + deg > target_vol) {
      break;
    }
    curvolume += deg;
    volume[i] = curvolume;
    if (curvolume == 0 || total_degree-curvolume==0) {
      conductance[i] = 1;
    } else {
      conductance[i] = (double)curcutsize/
                        (double)std::min(curvolume,total_degree-curvolume);
    }
    //printf("%5i : cut=%6i vol=%6i prval=%8g cond=%f\n", i, curcutsize, curvolume, it->second, conductance[i]);
  }
  // we stopped the iteration when it finished, or when it hit target_vol
  size_t lastind = i;
  double mincond = std::numeric_limits<double>::max();
  size_t mincondind = 0; // set to zero so that we only add one vertex 
  for (i=0; i<lastind; i++) {
    if (conductance[i] < mincond) {
      mincond = conductance[i];
      mincondind = i;
    }
  }
  //printf("mincond=%f mincondind=%i\n", mincond, mincondind);
  if (lastind == 0) {
    // add a case 
    mincond = 0.0;
  }
  i = 0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend && i<mincondind+1; ++it, ++i) {
    cluster.push_back(it->first);
  }
  if (outcond) { *outcond = mincond; }
  if (outvolume) { *outvolume = volume[mincondind]; }
}

struct local_pagerank_stats {
    double conductance;
    int volume;
    int support;
    int steps;
    double eps;
};
    

/** Cluster will contain a list of all the vertices in the cluster
 * @param start the vertex to use as the seed for the cluster
 * @param alpha the value of alpha in the PageRank computation
 * @param target_vol the approximate number of edges in the cluster
 * @param max_vol the maximum volume of a cluster (unused)
 * @param p the pagerank vector
 * @param r the residual vector
 * @param a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */
template <class Queue>
int hypercluster_pagerank_multiple(sparserow* G, 
    int start, double alpha, int target_vol, int max_vol,
    sparsevec& p, sparsevec &r, Queue& q,
    std::vector<int>& cluster, local_pagerank_stats *stats)
{
  // reset data
  p.map.clear();
  r.map.clear();
  q.empty();
  
  assert(target_vol > 0);
  assert(max_vol > 0);
  assert(alpha < 1.0); assert(alpha > 0.0);
  assert(start >= 0); assert(start < G->n); 
  
  r.map[start] = 1.0;
  //double pr_eps = 1.0/std::max((double)sr_degree(G,start)*(double)target_vol, 100.0);
  double pr_eps = 1.0/std::max(10.*(double)target_vol, 100.0);
  if (stats) { stats->eps = pr_eps; }
  
  //printf("find_cluster: start=%7i target_vol=%7i max_vol=%7i alpha=%5.3f pr_eps=%f\n", start, target_vol, max_vol, alpha, pr_eps);
  
  // calculate an integer number of maxsteps
  double maxsteps = 1./(pr_eps*(1.-alpha));
  maxsteps = std::min(maxsteps, 0.5*(double)std::numeric_limits<int>::max());
      
  int nsteps = compute_local_pagerank(G, r, p, alpha, pr_eps, (int)maxsteps, q);
  int support = r.map.size();
  if (stats) { stats->steps = nsteps; }
  if (stats) { stats->support = support; }

  // scale the probablities by their degree
  for (sparsevec::map_type::iterator it=p.map.begin(),itend=p.map.end();
    it!=itend;++it) {
    it->second *= 1.0/(double)std::max(sr_degree(G,it->first),1);  
  }
  double *outcond = NULL;
  int *outvolume = NULL;
  if (stats) { outcond = &stats->conductance; }
  if (stats) { outvolume = &stats->volume; }
  cluster_from_sweep(G, p, max_vol, cluster, outcond, outvolume);
  return (0);
}

int hypercluster_pagerank_single(sparserow* G, int start, double alpha, 
                                 int target_vol, int max_vol, 
                          std::vector<int>& cluster)
{
  sparsevec p, r;
  std::queue<int> q;
  return hypercluster_pagerank_multiple(G, 
          start, alpha, target_vol, max_vol, p, r, q, cluster, NULL);
}


int hypercluster_pagerank(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t overlap, size_t minsize,
                          double maxcond,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<double>& conductance,
                          std::vector<int>& order,
                          FILE* statsfile)
{
  std::vector<size_t> vertex_overlap(G->n);
  size_t nmarked1=0, nmarked2=0;
  sparsevec p, r;
  array_queue q(G->m);
  int nclusters = 0;
  double totalvol = 0;
  local_pagerank_stats stats;
  double t0=sf_time();
  
  for (int oi = 0; oi < G->m; oi++) {
    int i = order[oi];
    if (vertex_overlap[i]>=overlap) { continue; }
    int deg = sr_degree(G,i);
    if (deg == 0) { continue; }
    double curexpand = (double)(deg+1)*expand;
    // wow, be careful with the upper limit here, I had the 
    // random number generator generate a number less than 0.1 115 
    // times in a row, which caused target_vol to become negative.
    double maxexpand = (double)G->m;
    while (curexpand < maxexpand) { // while the loop hasn't broken...
      double rand_val = sf_rand();
      if (rand_val < 0.1) {
        break;
      }
      curexpand *= expandfactor;
      //printf("curexpand=%f rand_val=%f\n", curexpand, rand_val);
    }
    
    if (curexpand > maxexpand) { curexpand = maxexpand; }
    
    int target_vol = std::min(10*max_vol, (int)curexpand);
    
    assert(target_vol > 0);
    assert(target_vol <= G->m);
    
    std::vector<int> cluster;
    hypercluster_pagerank_multiple(G, i, 
      alpha, target_vol, max_vol, p, r, q,
      cluster, &stats);
    
    if (cluster.size()>=minsize && stats.conductance<=maxcond) {
      size_t nused = 0;
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        vertex_overlap[cluster[ci]]+=1;
        if (vertex_overlap[cluster[ci]] == overlap) {
          nmarked2++;
          nused++;
        }
        if (vertex_overlap[cluster[ci]] == 1) {
          nmarked1++;
          nused++;
        }
      }
      
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        clusters.push_back(std::make_pair(nclusters, cluster[ci]));
      }
      nclusters++;
      totalvol += stats.volume;
      
      conductance.push_back(1.0-stats.conductance);
      
      // code to double check that we get the correct conductance
      // double myconductance = (double)sr_graph_cutsize(G, &cluster[0], cluster.size())/
      //                           (double)sr_graph_volume(G, &cluster[0], cluster.size());
      // printf("mycond = %g\n", myconductance); 
      
      // vertices, volume, cond   steps, support   rate   utility, completion
      printf("%6zu %7i %5.3f   %6i %6i %8.2e  %7.1f  %5.1f%% %5.1f%%\n",
        cluster.size(), stats.volume, stats.conductance, 
        stats.steps, stats.support, stats.eps, totalvol/(sf_time()-t0),
        100.*(double)nused/(double)cluster.size(),
        100.*(double)nmarked2/(double)G->n);
        
      if (statsfile) {
        fprintf(statsfile, "%i %i %zu %i %g %i %i %i\n",
            i, sr_degree(G,i), cluster.size(), stats.volume, 
            stats.conductance, stats.steps, stats.support,
            target_vol);
      }
        
      //printf("nv: %8Zu vol: %8i cond: %6.4f (rate: %8.2f vol/sec) %4.1f%% %4.1f%% marked\n",
      //  cluster.size(), volume, cond, totalvol/(sf_time()-t0), 
      //  100.*(double)nused/(double)cluster.size(),
      //  100.*(double)nmarked2/(double)G->n);
    }
  }

  return (nclusters);
}


int hypercluster_pagerank_full(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t overlap, size_t minsize,
                          double maxcond,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<double>& conductance,
                          std::vector<int>& order,
                          FILE* statsfile)
{
  std::vector<size_t> vertex_overlap(G->n);
  size_t nmarked1=0, nmarked2=0;
  sparsevec p, r;
  array_queue q(G->m);
  int nclusters = 0;
  double totalvol = 0;
  local_pagerank_stats stats;
  double t0=sf_time();
  
  // start a small cluster from each vertex always
  printf("Initial clustering\n");
  for (int oi = 0; oi < G->m; oi++) {
    int i = order[oi];
    int deg = sr_degree(G,i);
    if (deg == 0) { continue; }
    double curexpand = (double)(deg+1);
    int target_vol = std::min(10*max_vol, (int)curexpand);
    
    assert(target_vol > 0);
    assert(target_vol <= G->m);
    
    std::vector<int> cluster;
    hypercluster_pagerank_multiple(G, i, 
      alpha, target_vol, max_vol, p, r, q,
      cluster, &stats);
    
    if (cluster.size()>=minsize && stats.conductance<=maxcond) {
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        clusters.push_back(std::make_pair(nclusters, cluster[ci]));
      }
      nclusters++;
      totalvol += stats.volume;
      
      conductance.push_back(1.0-stats.conductance);
      
      // vertices, volume, cond   steps, support   rate   utility, completion
      printf("%6zu %7i %5.3f   %6i %6i %8.2e  %7.1f  %5.1f%%\n",
        cluster.size(), stats.volume, stats.conductance, 
        stats.steps, stats.support, stats.eps, totalvol/(sf_time()-t0),
        100.*(double)nmarked2/(double)G->n);
        
      if (statsfile) {
        fprintf(statsfile, "%i %i %zu %i %g %i %i %i\n",
            i, sr_degree(G,i), cluster.size(), stats.volume, 
            stats.conductance, stats.steps, stats.support,
            target_vol);
      }
    }
  }
  
  // now try and find bigger clusters
  for (int oi = 0; oi < G->m; oi++) {
    int i = order[oi];
    if (vertex_overlap[i]>=overlap) { continue; }
    
    for (int tryi=0; tryi < overlap; ++tryi) {
      if (vertex_overlap[i]>=overlap) { continue; }
      int deg = sr_degree(G,i);
      if (deg == 0) { continue; }
      double curexpand = (double)(deg+1)*expand*pow(expandfactor,(double)(tryi+1));
      double maxexpand = (double)G->m;
      if (curexpand > maxexpand) { curexpand = maxexpand; }
      
      int target_vol = std::min(10*max_vol, (int)curexpand);
      
      assert(target_vol > 0);
      assert(target_vol <= G->m);
      
      std::vector<int> cluster;
      hypercluster_pagerank_multiple(G, i, 
        alpha, target_vol, max_vol, p, r, q,
        cluster, &stats);
      
      if (cluster.size()>=minsize && stats.conductance<=maxcond) {
        size_t nused = 0;
        for (size_t ci=0; ci<cluster.size(); ++ci) {
          vertex_overlap[cluster[ci]]+=1;
          if (vertex_overlap[cluster[ci]] == overlap) {
            nmarked2++;
            nused++;
          }
          if (vertex_overlap[cluster[ci]] == 1) {
            nmarked1++;
            nused++;
          }
        }
        
        for (size_t ci=0; ci<cluster.size(); ++ci) {
          clusters.push_back(std::make_pair(nclusters, cluster[ci]));
        }
        nclusters++;
        totalvol += stats.volume;
        
        conductance.push_back(1.0-stats.conductance);
        
        // vertices, volume, cond   steps, support   rate   utility, completion
        printf("%6zu %7i %5.3f   %6i %6i %8.2e  %7.1f  %5.1f%% %5.1f%%\n",
          cluster.size(), stats.volume, stats.conductance, 
          stats.steps, stats.support, stats.eps, totalvol/(sf_time()-t0),
          100.*(double)nused/(double)cluster.size(),
          100.*(double)nmarked2/(double)G->n);
          
        if (statsfile) {
          fprintf(statsfile, "%i %i %zu %i %g %i %i %i\n",
              i, sr_degree(G,i), cluster.size(), stats.volume, 
              stats.conductance, stats.steps, stats.support,
              target_vol);
        }
      }
    }
  }

  return (nclusters);
}



int multicluster_pagerank(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t minsize,
                          double maxcond,
                          size_t maxvisits,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<double>& conductance,
                          std::vector<int>& vertices,
                          FILE* statsfile)
{
  std::vector<size_t> vertex_visits;
  std::vector<size_t> vertex_overlap(G->n,0);
  size_t nmarked=0;
  sparsevec p, r;
  array_queue q(G->m);
  int nclusters = 0;
  double totalvol = 0;
  local_pagerank_stats stats;
  double t0=sf_time();
  
  if (maxvisits > 0) {
      vertex_visits.resize(G->n,0);
  }
  
  for (size_t si = 0; si < vertices.size(); ++si) {
    int i = vertices[si];
    
    if (maxvisits > 0 && vertex_visits[i] > maxvisits) {
        continue;
    }
    
    int deg = sr_degree(G,i);
    if (deg == 0) { continue; }
    double curexpand = (double)(deg+1)*expand;
    // wow, be careful with the upper limit here, I had the 
    // random number generator generate a number less than 0.1 115 
    // times in a row, which caused target_vol to become negative.
    double maxexpand = (double)G->m;
    while (curexpand < maxexpand) { // while the loop hasn't broken...
      double rand_val = sf_rand();
      if (rand_val < 0.1) {
        break;
      }
      curexpand *= expandfactor;
      //printf("curexpand=%f rand_val=%f\n", curexpand, rand_val);
    }
    
    if (curexpand > maxexpand) { curexpand = maxexpand; }
    
    int target_vol = std::min(10*max_vol, (int)curexpand);
    
    assert(target_vol > 0);
    assert(target_vol <= G->m);
    
    std::vector<int> cluster;
    hypercluster_pagerank_multiple(G, i, 
      alpha, target_vol, max_vol, p, r, q,
      cluster, &stats);
      
    if (maxvisits > 0) { // only update this if we are tracking visits
        for (sparsevec::map_type::const_iterator it = p.map.begin(); 
            it != p.map.end(); ++it) {
            vertex_visits[it->first] += 1;
        }
    }
    
    if (cluster.size()>=minsize && stats.conductance<=maxcond) {
      size_t nused = 0;
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        vertex_overlap[cluster[ci]]+=1;
        if (vertex_overlap[cluster[ci]] == 1) {
          nmarked++;
          nused++;
        }
      }
      
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        clusters.push_back(std::make_pair(nclusters, cluster[ci]));
      }
      nclusters++;
      totalvol += stats.volume;
      
      conductance.push_back(1.0-stats.conductance);
      
      // vertices, volume, cond   steps, support   rate   utility, completion
      printf("%6zu %7i %5.3f   %6i %6i %8.2e  %7.1f  %5.1f%% %5.1f%%\n",
        cluster.size(), stats.volume, stats.conductance, 
        stats.steps, stats.support, stats.eps, totalvol/(sf_time()-t0),
        100.*(double)nused/(double)cluster.size(),
        100.*(double)nmarked/(double)G->n);
        
      if (statsfile) {
        fprintf(statsfile, "%i %i %zu %i %g %i %i %i\n",
            i, sr_degree(G,i), cluster.size(), stats.volume, 
            stats.conductance, stats.steps, stats.support,
            target_vol);
      }
    
    }
  }

  return (nclusters);
}



/** 
 * This function performs repeated clusters of the graph as follows
 */
/*int hypercluster_pagerank_boundary(sparserow* G, double alpha, int max_vol,
                          double expand, double expandfactor,
                          size_t sweeps, size_t boundary,
                          std::vector<std::pair<int,int> >& clusters,
                          std::vector<int>& order)
{
  std::vector<size_t> vertex_overlap(G->n);
  size_t nmarked1=0, nmarked2=0;
  sparsevec p, r;
  array_queue q(G->m);
  int nclusters = 0;
  double totalvol = 0;
  local_pagerank_stats stats;
  double t0=sf_time();
  
  for (size_t sweep=0; sweep < sweeps; ++sweep) {
    std::vector<bool> picked(G->n, false);
    std::vector<bool> marked(G->n, false);
    
    // sweep 1: pick well separated clusters
    for (int oi = 0; oi < G->m; oi++) {
      int i = order[oi];
      if (marked[i]) { continue; }
      int deg = sr_degree(G,i);
      if (deg == 0) { continue; }
      double curexpand = (double)(deg+1)*expand;
      // wow, be careful with the upper limit here, I had the 
      // random number generator generate a number less than 0.1 115 
      // times in a row, which caused target_vol to become negative.
      double maxexpand = (double)G->m;
      while (curexpand < maxexpand) { // while the loop hasn't broken...
        double rand_val = sf_rand();
        if (rand_val < 0.1) {
          break;
        }
        curexpand *= expandfactor;
      }
      if (curexpand > maxexpand) { curexpand = maxexpand; }
      int target_vol = std::min(10*max_vol, (int)curexpand);
      assert(target_vol > 0);
      assert(target_vol <= G->m);
      std::vector<int> cluster;
      hypercluster_pagerank_multiple(G, i, 
        alpha, target_vol, max_vol, p, r, q,
        cluster, &stats);
  
  for (int oi = 0; oi < G->m; oi++) {
    int i = order[oi];
    if (vertex_overlap[i]>=overlap) { continue; }
    int deg = sr_degree(G,i);
    if (deg == 0) { continue; }
    double curexpand = (double)(deg+1)*expand;
    // wow, be careful with the upper limit here, I had the 
    // random number generator generate a number less than 0.1 115 
    // times in a row, which caused target_vol to become negative.
    double maxexpand = (double)G->m;
    while (curexpand < maxexpand) { // while the loop hasn't broken...
      double rand_val = sf_rand();
      if (rand_val < 0.1) {
        break;
      }
      curexpand *= expandfactor;
    }
    if (curexpand > maxexpand) { curexpand = maxexpand; }
    int target_vol = std::min(10*max_vol, (int)curexpand);
    assert(target_vol > 0);
    assert(target_vol <= G->m);
    std::vector<int> cluster;
    hypercluster_pagerank_multiple(G, i, 
      alpha, target_vol, max_vol, p, r, q,
      cluster, &stats);
    
    if (cluster.size()>1) {
      size_t nused = 0;
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        vertex_overlap[cluster[ci]]+=1;
        if (vertex_overlap[cluster[ci]] == overlap) {
          nmarked2++;
          nused++;
        }
        if (vertex_overlap[cluster[ci]] == 1) {
          nmarked1++;
          nused++;
        }
      }
      
      for (size_t ci=0; ci<cluster.size(); ++ci) {
        clusters.push_back(std::make_pair(nclusters, cluster[ci]));
      }
      nclusters++;
      totalvol += stats.volume;
      
      // vertices, volume, cond   steps, support   rate   utility, completion
      printf("%6zu %7i %5.3f   %6i %6i %8.2e  %7.1f  %5.1f%% %5.1f%%\n",
        cluster.size(), stats.volume, stats.conductance, 
        stats.steps, stats.support, stats.eps, totalvol/(sf_time()-t0),
        100.*(double)nused/(double)cluster.size(),
        100.*(double)nmarked2/(double)G->n);
        
      //printf("nv: %8Zu vol: %8i cond: %6.4f (rate: %8.2f vol/sec) %4.1f%% %4.1f%% marked\n",
      //  cluster.size(), volume, cond, totalvol/(sf_time()-t0), 
      //  100.*(double)nused/(double)cluster.size(),
      //  100.*(double)nmarked2/(double)G->n);
    }
  }

  return (nclusters);
}*/
