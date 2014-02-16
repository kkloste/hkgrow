/** @file hypercluster.cc
 * Routines to hypercluster a graph.  
 * A hyperclustering is a set of overlapping clusters.
 */

/*
 * David F. Gleich
 * Copyright, 2010.
 */


/** 
 * History
 * -------
 * :2008-09-01: Initial coding 
 * :2010-02-03: Fixed outputing isolated vertices by merging them with others
 * :2010-02-05: Switched to non-randomized by default
 * :2010-07-28: Added centrality option
 * :2010-09-29: Added larger overlap option
 * :2010-10-21: 
 * :2011-06-21: Added just individual runs
 */

#define NOMINMAX

#include <iostream>
#include <string>
#include <utility>
#include <limits>
#include <algorithm>
#include <queue>
#include <fstream>
#include "hypercluster.h"
#include "sparfun_util.h"
#include "sparvec.h"

#ifdef __APPLE__
#include <tr1/random>
#else
#include <random>
#endif

#pragma warning(disable:4996)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <assert.h>

/** Convert a path into a basename and an extension
 * The extension of a path is the string following
 * the last period in the path.
 *
 * @example 
 * path_basename_and_extension("../data/myfile.ext") 
 *   returns std::pair("../data/myfile","ext")
 * path_basename_and_extension("../data/myfile.label.txt") 
 *   returns std::pair("../data/myfile.label","txt")
 * path_basename_and_extension("../data/myfile") 
 *   returns std::pair("../data/myfile.label","")
 *
 * @param p the path
 * @return a pair of (basename, extension)
 */
std::pair<std::string,std::string> path_basename_and_extension(const char* p)
{
    std::string s = p;
    std::string::size_type lastdot = s.rfind(".");
    if (lastdot == std::string::npos) {
        return make_pair(s, std::string(""));
    } else {
        return make_pair(s.substr(0, lastdot), s.substr(lastdot+1));
    }
}

void usage(FILE* s, bool exitprog=false) {
  const char* usagestr = 
    "hypercluster.exe graphfile [-out clusterfilename]\n"
    "  [-maxvol <int>|<float>] [-method <methodname>] [-seed <int>]\n"
    "  [additional method options, get a list with --help]\n"
    "\n"
    "hypercluster tries to produce clusters with useful overlap.\n"
    "\n"
    " -out clusterfilename : defaults to graphname.hcluster.<maxvol>\n"
    " -maxvol <int> : uses the positive integer (>1) as the maximum cluster volume\n"
    " -maxvol <float> : uses the floating value [0,1] to set the maximun\n"
    "                   cluster volume to <float>*number of edges\n"
    " -seed <int> : control randomness with a seed, -1 => use time\n"
    " -method <methodname> : changes the method to produce the clusters\n"
    " -two_core : only operate on the two core and reattach tree-like\n"
    "             pieces at the end\n"
    "\n"
    " -centrality <centrality_filename> : a file with centrality scores for all nodes\n"
    "\n"
    "Defaults:\n"
    "  maxvol = max(0.01*edges,100), method = pagerank, seed = 0 (non-random)\n"
    "  centrality = <None>\n"
    ;
  
  fprintf(s,"%s\n",usagestr);
  if (exitprog) { exit(1); }
}

int default_max_vol(int nedges) {
  return std::max((int)(0.01*nedges),100);
}

void hypercluster_pagerank_usage(FILE *s, bool exitprog=false) {
  const char* usagestr = 
    "pagerank hyperclustering options\n"
    " -alpha <float> : the value of alpha in PageRank (0,1)\n"
    " -expand <int> : an expansion factor for a cluster (>=1)\n"
    " -expand <float> : a starting expansion factor for a cluster\n"
    "                   (of total edges) (>1)\n"
    " -expandfactor <float> : a random expansion factor for a cluster (>1)\n";
    " -overlap <int> : the desired overlap on each vertex\n"
    "\n"
    "Defaults:\n"
    "  alpha = 0.99; expand = 3; expandfactor = 1.1; overlap = 2\n"
    ;
  fprintf(s,"%s\n",usagestr);
  if (exitprog) { exit(1); }
}

void hypercluster_pagerank_indep_usage(FILE *s, bool exitprog=false) {
  const char* usagestr = 
    "pagerank hyperclustering options\n"
    " REQUIRED\n"
    " -vertices <file>: a list of valid vertex identifers to start clusters from\n"
    " OPTIONS\n"
    " -alpha <float> : the value of alpha in PageRank (0,1)\n"
    " -expand <int> : an expansion factor for a cluster (>=1)\n"
    " -expand <float> : a starting expansion factor for a cluster\n"
    "                   (of total edges) (>1)\n"
    " -expandfactor <float> : a random expansion factor for a cluster (>1)\n"
    
    " -maxcond <float> : a value for the maximum conductance to output\n"
    "\n"
    "Defaults:\n"
    "  alpha = 0.99; expand = 3; expandfactor = 1.1;\n"
    ;
  fprintf(s,"%s\n",usagestr);
  if (exitprog) { exit(1); }
}

bool file_write_test(const char* filename) {
  FILE *f = fopen(filename, "a");
  bool rval = false;
  if (f) {
    return true;
  } 
  fclose(f);
  return rval;
}

bool read_centrality(const char* filename, 
        std::vector<float>& centrality, int nnodes) {
  int nv = 0;
  std::ifstream f(filename);
  while (f) {
    f >> centrality[nv];
    nv++;
    if (nv>=nnodes) {
      break;
    }
  }
  if (nv != nnodes) {
    fprintf(stderr,"error reading centrality file %s : only %i values of %i read\n", 
      filename, nv, nnodes);
    return false;
  }
  return true;
}

bool read_vertex_list(const char* filename, 
        std::vector<int>& list, int nnodes) {
  int nv = 0;
  std::ifstream f(filename);
  while (f) {
    int v;
    f >> v;
    if (!f) {
        break;
    }
    list.push_back(v);
    nv++;
    if (nv==nnodes+1) {
      fprintf(stderr,"warning, reading many vertices from vertex list\n");
    }
    assert(v >= 0 && v < nnodes);
  }
  return true;
}


template <typename T, typename index_type = int>
class sort_order_comparison {
  const std::vector<T>& items;
public:
  sort_order_comparison(const std::vector<T>& i) : items(i) {}
  bool operator() (index_type i, index_type j) {
    return items[i] > items[j];
  }
};

/**
 * This function sorts in descending order
 * @param order this parameter is an output and will be initialized
 */
template <typename T>
void sort_permutation(const std::vector<T>& a, std::vector<int>& order) {
  order.resize(a.size());
  for (size_t i=0; i<a.size(); ++i) {
    order[i] = (int)i;
  }
  sort_order_comparison<T> comp(a);
  std::sort( order.begin(), order.end(), comp );
}

template <typename Index>
void inverse_permutation(const std::vector<Index>& perm, std::vector<Index>& inv) {
  inv.resize(perm.size());
  for (size_t i=0; i<perm.size(); ++i) {
    inv[perm[i]] = i;
  }
}

/** Check for isolated vertices (not in any clusters)
 * @param n the size of the grpah
 * @param clusters a list of cluster,vertex pairs
 */
bool check_isolated_vertices( int n, 
    std::vector< std::pair<int, int> >& clusters,
    std::vector<int>& flagged) 
{
  for (int i=0; i<n; i++) { flagged[i] = 0; }
  int nflagged = 0;
  typedef std::vector< std::pair<int, int> >::const_iterator iter;
  for (iter i=clusters.begin(); i!=clusters.end(); i++) {
    int vertex = i->second;
    if (!flagged[vertex]) {
      flagged[vertex]=1;
      nflagged ++;
    }
  }
  if (nflagged != n) {
    return true;
  } else {
    return false;
  }
}

/** Build a map from vertices to clusters to help fix isolated vertices
 */
sparserow* index_clusters(std::vector< std::pair<int, int> >& clusters,
    int nclusters, int nverts) 
{    
  typedef std::vector< std::pair<int, int> >::const_iterator iter;
  sparserow* s = sparserow_fullalloc(nverts, clusters.size(), false);
  if (s) {
    s->m = nverts; 
    s->n = nclusters;
    memset(s->ai, 0, sizeof(int)*(s->m+1));
    int i, j, nzi;
    for (iter ci=clusters.begin(); ci!=clusters.end(); ++ci) {
      s->ai[ci->second+1]++; 
    }
    for (nzi=0, j=0; j<s->m+1; j++) {
      s->ai[j] = (nzi += s->ai[j]);
    }
    for (iter ci=clusters.begin(); ci!=clusters.end(); ++ci) {
      s->aj[s->ai[ci->second]] = ci->first;
      s->ai[ci->second]++; 
    }
    for (j=s->m; j>0; j--) {
      s->ai[j] = s->ai[j-1];
    }
    s->ai[0] = 0;
  }
  return s;
}  
  
/** Make a map from two_core vertices to one_core vertices that touch them.
 * 
 * Suppose that vertices 1 and 2 are in a two-node chain from vertex 3.
 * e.g. 3-2-1
 * but that vertex 3 has many other connections.   
 * Then the output from this function will assign 2 and 1 to vertex 3.
 *  child_list[3] = [2,1]
 * 
 * @param r the graph
 * @param core_map the core numbers of each vertex
 * @param child_list an array of arrays to store the 
 *   one-core vertices for any vertex in the two core
 */
void build_one_core_child_list(const sparserow* r, int *core_map, 
    std::vector< std::vector<int> >& child_list) 
{
  // set each vertex to unassigned initially
  std::vector<int> assigned(r->n, -1);
  
  
  std::vector<int> path; // the list of one-core vertices on our path
  std::queue<int> queue; // the queue for the one-core bfs
  for (int i=0; i<r->n; ++i) {
    // we only start the bfs from vertices in the two-core
    if (core_map[i]<=1) {
      continue;
    }
    // from this vertex, we crawl for all vertices in the one-core.
    path.clear();
    assert(queue.size() == 0);
    for (int nzi=r->ai[i]; nzi<r->ai[i+1]; ++nzi) {
      int j = r->aj[nzi];
      if (core_map[j] <= 1) {
        queue.push(j);
      }
    }
    
    while (queue.size() > 0) {
      int vertex = queue.front();
      queue.pop();
      assert(core_map[vertex] <= 1);
      path.push_back(vertex);
      assigned[vertex] = i;
      
      for (int nzi=r->ai[vertex]; nzi<r->ai[vertex+1]; ++nzi) {
        int j = r->aj[nzi];
        if (core_map[j] <= 1) {
          if (assigned[j] == -1) {
            queue.push(j);
          } else {
            assert(assigned[j] == i); // all other one core must be assigned to us too
          }
        } else {
          assert(j == i); // the only two core vertex we should see is i
        }
      }
    }
    
    // save the child list
    child_list[i] = path;
  }
}    

/** Check for isolated vertices in a clustering and fix them if found
 * This routine first checks if any vertices have NOT been assigned to a
 * cluster.  If so, then it builds a map from vertices to clusters.
 * Using this map, it first tries to map the vertex into the 
 * most common neighbor.
 * @param r the graph
 * @param nisolated the number of isolated vertices
 * @return the number of new clusters
 */
int check_and_fix_isolated_vertices( sparserow* r, 
    std::vector< std::pair<int, int> >& clusters, int nclusters, 
    size_t *nisolated) 
{
  std::vector<int> isolated_flag(r->n);
  if (nisolated) { *nisolated = 0; }
  if (check_isolated_vertices(r->n, clusters, isolated_flag) == false) {
    return 0; // no isolated vertices
  }
  // fix the isolated vertices
  // map vertices to clusters, arg... annoying that we have to build this here
  sparserow* c = index_clusters(clusters, nclusters, r->n);
  
  bool newcluster = false;
  
  // look over all vertices to see the connectivity
  for (size_t vi=0; vi<r->n; ++vi) {
    // is vertex i isolated
    if (isolated_flag[vi] == 0) {
      if (nisolated) { (*nisolated) ++; }  
      if (sr_degree(r, vi)==0) {
        // assign to any cluster, but lets make a new one of isolated vertices
        newcluster = true;
        isolated_flag[vi] = 1;
        clusters.push_back(std::make_pair(nclusters, vi));
      } else {
        // the graph is symmetric, so check the concensus vote
        // of the clusters of the neighbors
        sparsevec neighbor_clusters;
        for (int nzi=r->ai[vi]; nzi<r->ai[vi+1]; ++nzi) {
          int ni=r->aj[nzi];
          // look at all clusters the neighbor is in
          for (int nzj=c->ai[ni]; nzj<c->ai[ni+1]; ++nzj) {
            int ci = c->aj[nzj]; // the cluster index
            neighbor_clusters.map[ci] = neighbor_clusters.get(ci)+1.;
          }
        }
        if (neighbor_clusters.map.size() == 0) {
          // none of the neighbors were in any cluster, so just
          // assign them to our new cluster
          newcluster = true;
          clusters.push_back(std::make_pair(nclusters, vi));
        } else {
          // find the most frequent element
          int ci = neighbor_clusters.max_index();
          // assign the vertex to this cluster
          clusters.push_back(std::make_pair(ci, vi));
        }
      }
    }
  }
  // return the number of new clusters added.
  int rval;
  if (newcluster) { 
    rval = 1;
  } else {
    rval = 0;
  }
  free_sparse(c); free(c);
  return rval;
}  
  
        

/**
 * @param maxvol = -0 => use default_max_vol()
 */
int hypercluster_pagerank(const char* graphfilename, double maxvol, 
                          const char* outputfilebase,
                          int argc, char **argv)
{
  double alpha = 0.99;
  double expand = 3;
  double expandfactor = 1.1;
  const char* centrality_filename = NULL;
  int overlap = 2;
  bool two_core = false;
  bool full = false;
  double maxcond = 0.55;
  size_t minsize = 10;
  const char* stats_filename = NULL;

  // parse our options
  for (int argi=2; argi < argc; ++argi) {
    if (strcmp(argv[argi],"-alpha")==0) {
      if (argi+1<argc) {
        double alphaval = atof(argv[++argi]);
        if (alphaval >= 1 || alphaval <= 0) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-alpha", alphaval, "outside range (0,1)");
          return (1);
        } else {
          alpha = alphaval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n", 
          "-alpha", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-expand")==0) {
      if (argi+1<argc) {
        double expandval = atof(argv[++argi]);
        if (expandval <= 0) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-expand", expandval, "outside range (0,inf)");
          return (1);
        } else {
          expand = expandval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n", 
          "-expand", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-centrality")==0) {
      if (argi+1<argc) {
        centrality_filename = argv[++argi];
      } else {
        fprintf(stderr,"error parsing argument %s : %s\n", 
          "-centrality", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-statsfile")==0) {
      if (argi+1<argc) {
        stats_filename = argv[++argi];
      } else {
        fprintf(stderr,"error parsing argument %s : %s\n", 
          "-statsfile", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-expandfactor")==0) {
      if (argi+1<argc) {
        double expandfactorval = atof(argv[++argi]);
        if (expandfactorval < 1) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-expandfactor", expandfactorval, "outside range [1,inf)");
          return (1);
        } else {
          expandfactor = expandfactorval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n", 
          "-expandfactor", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-overlap")==0) {
      if (argi+1<argc) {
        int overlapval = atoi(argv[++argi]);
        if (overlapval < 1) {
          fprintf(stderr,"error parsing pagerank argument %s : value %i %s\n",
            "-overlap", overlapval, "outside range (1,inf)");
          return (1);
        } else {
          overlap = overlapval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-overlap", "no value supplied");
      }
    } else if (strcmp(argv[argi],"-minsize")==0) {
      if (argi+1<argc) {
        int minsizeval = atoi(argv[++argi]);
        if (minsizeval < 1) {
          fprintf(stderr,"error parsing pagerank argument %s : value %i %s\n",
            "-minsize", minsizeval, "outside range (1,inf)");
          return (1);
        } else {
          minsize = (size_t)minsizeval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-minsize", "no value supplied");
      }
    } else if (strcmp(argv[argi],"-maxcond")==0) {
      if (argi+1<argc) {
        double maxcondval = atof(argv[++argi]);
        if (maxcondval > 1. || maxcondval < 0.) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-maxcond", maxcondval, "outside range (0,1)");
          return (1);
        } else {
          maxcond = (double)maxcondval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-maxcond", "no value supplied");
      }
    } else if (strcmp(argv[argi],"-two_core")==0) {
      two_core = true;
    } else if (strcmp(argv[argi],"-full")==0) {
      full = true;
    } else if (strcmp(argv[argi],"--help")==0) {
      usage(stdout,true);
    } else if (strlen(argv[argi])>0) {
      fprintf(stderr,"unknown pagerank option %s\n", argv[argi]);
      //return (1);
    }
  }

  using namespace std;

  // read the graph
  sparserow As = {0};
  sparserow *A = &As;

  double st = 0.0;  

  st = sf_time();
  int rval = load_generic(argv[1], A, false);

  if (rval != 0) {
    cout << "*** Error: " << rval << endl;
    return (rval); 
  }

  cout << "Loaded " << argv[1] << " in " << sf_time() - st << " secs" << endl;
  cout << "    nodes = " << A->m << endl;
  cout << "    edges = " << sr_nedges(A) << endl;
  cout << endl;

  st = sf_time();
  sparserow *G = sf_sym_or(A);
  if (G==NULL || sf_pack(G)!=0) {
    cout << "*** Error: symmetrizing" << endl;
    return (-1);
  } else {
    free_sparse(A);
  }

  cout << "Packed " << argv[1] << " in " << sf_time() - st << " secs" << endl;
  cout << "    nodes = " << G->m << endl;
  cout << "    edges = " << sr_nedges(G) << endl;
  cout << endl;
  
  
  
  // load centrality  and compute order
  std::vector<int> order;
  {
    std::vector<float> centrality(G->m);
    if (centrality_filename == NULL) {
      for (int i=0; i<G->m; ++i) {
        centrality[i] = (float)sr_degree(G,i);
      }
    } else {
      if (!read_centrality(centrality_filename, centrality, G->m)) {
        return (1);
      }
    }
    // build the order for centrality
    sort_permutation(centrality, order);
  }

  cout << "Parameters" << endl;
  int max_vol = 0;
  if (maxvol == -0) {
    max_vol = default_max_vol(sr_nedges(G));
  } else if (maxvol <= 1) {
    max_vol = (int)(maxvol*sr_nedges(G));
  } else {
    max_vol = (int)maxvol;
  }
  if (expand < 1) {
    expand = expand * (double)sr_nedges(G);
  }
  cout << "        maxvol = " << max_vol << endl;
  cout << "        expand = " << expand << endl;
  cout << "  expandfactor = " << expandfactor << endl;
  cout << "         alpha = " << alpha << endl;
  cout << "       overlap = " << overlap << endl;
  cout << "      two_core = " << two_core << endl;
  cout << "       minsize = " << minsize << endl;
  cout << "       maxcond = " << maxcond << endl;
  if (centrality_filename == NULL) {
    cout << "    centrality = " << "<DEGREE>" << endl;
  } else {
    cout << "    centrality = " << centrality_filename << endl;
  }
  if (stats_filename != NULL) {
    cout << "     statsfile = " << stats_filename << endl;
  }
  cout << endl;


  if (max_vol < 0) {
    cout << "Error: maxvol = " << maxvol << endl;
    cout << "It should be bigger than 0!" << endl;
    return (1);
  } else if (max_vol < 10) {
    cout << "    Warning: maxvol = " << maxvol << " is very small" << endl;
  }
  


  // test our output here to avoid doing computation only to have it fail
  std::string outfilename = outputfilebase;
  if (!file_write_test(outfilename.c_str())) {
    cout << endl;
    cout << " ** Error: couldn't open " << outfilename << endl;
    cout << endl;
    return 0;
  }
  
  FILE *statsfile = NULL;
  if (stats_filename) {
    statsfile = fopen(stats_filename, "wt");
    if (!statsfile) {
      cout << endl;
      cout << " ** Error: couldn't open " << stats_filename << endl;
      cout << endl;
      return 0;
    }
  }
  
  std::vector< std::pair<int, int> > clusters;
  std::vector< double > cscore;
  int nclusters;
  
  if (two_core) {
    assert(full == false);
    cout << "Using two-core" << std::endl;
    std::vector<int> core_numbers(G->n);
    sr_graph_core_numbers(G, &core_numbers[0]);
    
    // build two_core_parent_map
    // for all vertices if the two_core, this is the vertex
    // id itself, for all vertices in the one_core, it is the 
    // parent of the vertex that is in the two_core
    std::vector< std::vector<int> > child_list(G->n);
    build_one_core_child_list(G, &core_numbers[0], child_list);
    
    
    
    // extract the two-core-graph
    std::vector<size_t> set;
    std::vector<int> filter(G->n);
    size_t nset=0;
    
    for (int i=0; i<G->n; ++i) {
        if (core_numbers[i] > 1) {
            filter[i]=nset;
            set.push_back(i);
            nset+=1;
        } else {
            filter[i]=-1;
        }
    }
    assert(nset == set.size());
    sparserow *G2 = sf_rowcol_subset(G, &set[0], nset);
    // we can't free G because it's used to fix isolated vertices later
    
    cout << "  two-core vertices: " << nset << endl;
    // fix the order now too
    std::vector<int> order2;
    order2.reserve(nset);
    for (size_t i=0; i<order.size(); i++) {
        if (filter[order[i]]>=0) {
            order2.push_back(filter[order[i]]);
        }
    }
    
    // run the clustering
    nclusters = hypercluster_pagerank(G2, 
      alpha, max_vol, expand, expandfactor, 
      overlap, minsize, maxcond, 
      clusters, cscore, order2, statsfile);
      
    // now fix-up the clustering
    // we store last_assignment, because we will append things
    // to the list of clusters that we don't want to change
    // later
    size_t last_assignment = clusters.size();
    for (size_t i=0; i<last_assignment; i++) {
      clusters[i].second = set[clusters[i].second];
      int vertex = clusters[i].second;
      if (child_list[vertex].size() > 0) {
        size_t nadded = 0;
        // output these additional assignments
        for (size_t j=0; j<child_list[vertex].size(); ++j) {
          clusters.push_back(std::make_pair( clusters[i].first,
            child_list[vertex][j] ));
          nadded += 1;
        }
        //printf("twocore: added %Zi vertices to cluster %i\n", nadded, clusters[i].first);
      }
    }

  } else {

    if (full) {
      nclusters = hypercluster_pagerank_full(G, 
        alpha, max_vol, expand, expandfactor, 
        overlap, minsize, maxcond, 
        clusters, cscore, order, statsfile);
    } else {
      nclusters = hypercluster_pagerank(G, 
        alpha, max_vol, expand, expandfactor, 
        overlap, minsize, maxcond, 
        clusters, cscore, order, statsfile);
    }

  }
  
  if (nclusters <= 0) {
    cout << "    Warning: no clusters generated" << endl;
  }
  
  // sort clusters by cscore (which is really 1.0 - conductance)
  // this function sorts in descending order, and so gives us 
  // clusters by increasing conductance
  std::vector<int> corder;
  sort_permutation(cscore, corder);
  std::vector<int> cperm;
  inverse_permutation(corder, cperm);
  // now apply the permutation
  {
    typedef std::vector< std::pair<int, int> >::iterator it_type;
    for (it_type it=clusters.begin(),itend=clusters.end();it!=itend;++it) {
      it->first = cperm[it->first];
    }
  }


  size_t nisolated=0;
  nclusters += check_and_fix_isolated_vertices(
                            G, clusters, nclusters, &nisolated);
  if (nisolated) {
    cout << "  Fixed " << nisolated << " isolated vertices " << endl;
  }
   
  FILE *outfile = fopen(outfilename.c_str(),"wt");
  if (outfile) {
    fprintf(outfile,"%u %u %zu\n",nclusters,G->m,clusters.size());
    typedef std::vector< std::pair<int, int> >::const_iterator it_type;
    for (it_type it=clusters.begin(),itend=clusters.end();it!=itend;++it) {
      fprintf(outfile,"%u %u 1\n",it->first,it->second);
    }
  } else {
    cerr << "*** Error: Couldn't write output to " << outfilename << endl;
    return (-1);
  }

  return (0);
}


/**
 * @param maxvol = -0 => use default_max_vol()
 */
int multicluster_pagerank_indiv(const char* graphfilename, double maxvol, 
                          const char* outputfilebase,
                          int argc, char **argv)
{
  double alpha = 0.99;
  double expand = 3;
  double expandfactor = 1.1;
  const char* vertex_filename = NULL;
  double maxcond = 0.55;
  size_t minsize = 10;
  //int nsets = 5000;
  size_t maxvisits = 0;
  const char* stats_filename = NULL;

  // parse our options
  for (int argi=2; argi < argc; ++argi) {
    if (strcmp(argv[argi],"-alpha")==0) {
      if (argi+1<argc) {
        double alphaval = atof(argv[++argi]);
        if (alphaval >= 1 || alphaval <= 0) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-alpha", alphaval, "outside range (0,1)");
          return (1);
        } else {
          alpha = alphaval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n", 
          "-alpha", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-expand")==0) {
      if (argi+1<argc) {
        double expandval = atof(argv[++argi]);
        if (expandval <= 0) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-expand", expandval, "outside range (0,inf)");
          return (1);
        } else {
          expand = expandval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n", 
          "-expand", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-vertices")==0) {
      if (argi+1<argc) {
        vertex_filename = argv[++argi];
      } else {
        fprintf(stderr,"error parsing argument %s : %s\n", 
          "-vertices", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-statsfile")==0) {
      if (argi+1<argc) {
        stats_filename = argv[++argi];
      } else {
        fprintf(stderr,"error parsing argument %s : %s\n", 
          "-statsfile", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-expandfactor")==0) {
      if (argi+1<argc) {
        double expandfactorval = atof(argv[++argi]);
        if (expandfactorval < 1) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-expandfactor", expandfactorval, "outside range [1,inf)");
          return (1);
        } else {
          expandfactor = expandfactorval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n", 
          "-expandfactor", "no value supplied");
        return (1);
      }
    } else if (strcmp(argv[argi],"-minsize")==0) {
      if (argi+1<argc) {
        int minsizeval = atoi(argv[++argi]);
        if (minsizeval < 1) {
          fprintf(stderr,"error parsing pagerank argument %s : value %i %s\n",
            "-minsize", minsizeval, "outside range (1,inf)");
          return (1);
        } else {
          minsize = (size_t)minsizeval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-minsize", "no value supplied");
      }
    } else if (strcmp(argv[argi],"-maxvisits")==0) {
      if (argi+1<argc) {
        int maxvisitsval = atoi(argv[++argi]);
        if (maxvisitsval < 0) {
          fprintf(stderr,"error parsing pagerank argument %s : value %i %s\n",
            "-maxvisits", maxvisitsval, "outside range (0,inf)");
          return (1);
        } else {
          maxvisits = (size_t)maxvisitsval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-minsize", "no value supplied");
      }
    } else if (strcmp(argv[argi],"-maxcond")==0) {
      if (argi+1<argc) {
        double maxcondval = atof(argv[++argi]);
        if (maxcondval > 1. || maxcondval < 0.) {
          fprintf(stderr,"error parsing pagerank argument %s : value %f %s\n",
            "-maxcond", maxcondval, "outside range (0,1)");
          return (1);
        } else {
          maxcond = (double)maxcondval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-maxcond", "no value supplied");
      }
    } 
    /*else if (strcmp(argv[argi],"-nsets")==0) {
      if (argi+1<argc) {
        int nsetsval = atoi(argv[++argi]);
        if (nsetsval < 1) {
          fprintf(stderr,"error parsing pagerank argument %s : value %i %s\n",
            "-nsets", nsetsval, "outside range (1,inf)");
          return (1);
        } else {
          nsets = nsetsval;
        }
      } else {
        fprintf(stderr,"error parsing pagerank argument %s : %s\n",
          "-nsets", "no value supplied");
      }
    }*/ else if (strcmp(argv[argi],"--help")==0) {
      usage(stdout,true);
    } else if (strlen(argv[argi])>0) {
      fprintf(stderr,"unknown pagerank option %s\n", argv[argi]);
      //return (1);
    }
  }

  using namespace std;

  // read the graph
  sparserow As = {0};
  sparserow *A = &As;

  double st = 0.0;  

  st = sf_time();
  int rval = load_generic(argv[1], A, false);

  if (rval != 0) {
    cout << "*** Error: " << rval << endl;
    return (rval); 
  }

  cout << "Loaded " << argv[1] << " in " << sf_time() - st << " secs" << endl;
  cout << "    nodes = " << A->m << endl;
  cout << "    edges = " << sr_nedges(A) << endl;
  cout << endl;

  st = sf_time();
  sparserow *G = sf_sym_or(A);
  if (G==NULL || sf_pack(G)!=0) {
    cout << "*** Error: symmetrizing" << endl;
    return (-1);
  } else {
    free_sparse(A);
  }

  cout << "Packed " << argv[1] << " in " << sf_time() - st << " secs" << endl;
  cout << "    nodes = " << G->m << endl;
  cout << "    edges = " << sr_nedges(G) << endl;
  cout << endl;
  
  
  
  // load vertices
  std::vector<int> vertices;
  {
    if (vertex_filename == NULL) {
      for (int i=0; i<G->m; ++i) {
        vertices.push_back(i);
      }
    } else {
      if (!read_vertex_list(vertex_filename, vertices, G->m)) {
        return (1);
      }
    }
  }
  
  cout << "Read " << vertices.size() << " vertex ids" << endl;
  

  cout << "Parameters" << endl;
  int max_vol = 0;
  if (maxvol == -0) {
    max_vol = default_max_vol(sr_nedges(G));
  } else if (maxvol <= 1) {
    max_vol = (int)(maxvol*sr_nedges(G));
  } else {
    max_vol = (int)maxvol;
  }
  if (expand < 1) {
    expand = expand * (double)sr_nedges(G);
  }
  cout << "        maxvol = " << max_vol << endl;
  cout << "        expand = " << expand << endl;
  cout << "  expandfactor = " << expandfactor << endl;
  cout << "         alpha = " << alpha << endl;
  cout << "       minsize = " << minsize << endl;
  cout << "     maxvisits = " << maxvisits << endl;
  cout << "       maxcond = " << maxcond << endl;
  //cout << "         nsets = " << nsets << endl;
  if (vertex_filename == NULL) {
    cout << "      vertices = " << "<ALL>" << endl;
  } else {
    cout << "      vertices = " << vertex_filename << endl;
  }
  if (stats_filename != NULL) {
    cout << "     statsfile = " << stats_filename << endl;
  }
  cout << endl;


  if (max_vol < 0) {
    cout << "Error: maxvol = " << maxvol << endl;
    cout << "It should be bigger than 0!" << endl;
    return (1);
  } else if (max_vol < 10) {
    cout << "    Warning: maxvol = " << maxvol << " is very small" << endl;
  }
  
  // test our output here to avoid doing computation only to have it fail
  std::string outfilename = outputfilebase;
  if (!file_write_test(outfilename.c_str())) {
    cout << endl;
    cout << " ** Error: couldn't open " << outfilename << endl;
    cout << endl;
    return 0;
  }
  
  FILE *statsfile = NULL;
  if (stats_filename) {
    statsfile = fopen(stats_filename, "wt");
    if (!statsfile) {
      cout << endl;
      cout << " ** Error: couldn't open " << stats_filename << endl;
      cout << endl;
      return 0;
    }
  }
  
  std::vector< std::pair<int, int> > clusters;
  std::vector< double > cscore;
  int nclusters;
  
  
    
  nclusters = multicluster_pagerank(G, 
    alpha, max_vol, expand, expandfactor, 
    minsize, maxcond, maxvisits,
    clusters, cscore, vertices, statsfile);

  if (nclusters <= 0) {
    cout << "    Warning: no clusters generated" << endl;
  }
   
  FILE *outfile = fopen(outfilename.c_str(),"wt");
  if (outfile) {
    fprintf(outfile,"%u %u %zu\n",nclusters,G->m,clusters.size());
    typedef std::vector< std::pair<int, int> >::const_iterator it_type;
    for (it_type it=clusters.begin(),itend=clusters.end();it!=itend;++it) {
      fprintf(outfile,"%u %u 1\n",it->first,it->second);
    }
  } else {
    cerr << "*** Error: Couldn't write output to " << outfilename << endl;
    return (-1);
  }

  return (0);
}

int parse_standard_options_error(const char* arg, const char* message) {
  fprintf(stderr,"error parsing command line option %s: %s\n",arg,message);
  exit(1);
}
int parse_standard_options(int argc, char **argv,
    std::string& clusterfilename, double& maxvol, unsigned long& seed,
    std::string& methodname) {
  
  clusterfilename = ""; 
  maxvol = -0;
  //seed = 0; 
  // note that seed is initally set to time
  // and this function shouldn't change it unless
  // asked to by the user
  methodname = "pagerank";

  for (int argi=2; argi<argc; ++argi) {
    if (strcmp(argv[argi],"-out")==0) {
      if (argi+1 < argc) {
        clusterfilename = argv[++argi];
        argv[argi-1][0] = '\0';
        argv[argi][0] = '\0';
      } else {parse_standard_options_error("-out","missing output filename");}
    } else if (strcmp(argv[argi],"-maxvol")==0) {
      if (argi+1 < argc) {
        maxvol = atof(argv[argi+1]);
        argi++;
        argv[argi-1][0] = '\0';
        argv[argi][0] = '\0';
        if (errno == ERANGE) {
          parse_standard_options_error("-maxvol","cannot parse maxvol value");
        } else if (maxvol < 0) {
          parse_standard_options_error("-maxvol","maxvol must be >= 0");
        } else if (floor(maxvol) != maxvol && maxvol > 1) {
          printf("Warning: rounding maxval=%f to %i\n", 
            maxvol,(int)floor(maxvol));
          maxvol = floor(maxvol);
        }
      } else {parse_standard_options_error("-maxvol","missing volume");}
    } else if (strcmp(argv[argi],"-method") == 0) {
      if (argi+1 < argc) {
        methodname = argv[++argi];
        argv[argi-1][0] = '\0';
        argv[argi][0] = '\0';
      } else {parse_standard_options_error("-method","missing method name");}
    } else if (strcmp(argv[argi],"-seed") == 0) {
      if (argi+1 < argc) {
#ifdef __MSVC__          
        __int64 seed64 = _atoi64(argv[++argi]);
#else
        long long seed64 = atoll(argv[++argi]);
#endif 
        std::cout << "seed64 = " << seed64 << std::endl;        
        if (seed64<-1||seed64>std::numeric_limits<unsigned long>::max()) {
          fprintf(stderr, "seed %lli is outside the valid range (-1,%lu)\n",
                  seed64, std::numeric_limits<unsigned long>::max());
        }
        if (seed64 != -1) { seed = (unsigned long)seed64; }
        argv[argi-1][0] = '\0';
        argv[argi][0] = '\0';
      } else {parse_standard_options_error("-seed","missing seed value");}
    }
  }
  return (0);
}

int main(int argc, char **argv) 
{
  if (argc < 2) {
    usage(stdout, true); 
  }

  const char* graphfilename = argv[1];
  std::string clusterfilename;
  std::string methodname;

  double maxvol;
  unsigned long seed = (unsigned long)time(NULL);

  parse_standard_options(argc, argv, clusterfilename, 
    maxvol, seed, methodname);


  std::pair<std::string,std::string> graphname_parts = 
    path_basename_and_extension(graphfilename);

  std::cout << "Using seed = " << seed << std::endl;
  sf_srand(seed);

  std::string basefilename = graphname_parts.first;

  if (clusterfilename.empty()) {
    clusterfilename = basefilename + ".hcluster";
  }

  if (methodname == "pagerank") {
    return hypercluster_pagerank(graphfilename, 
      maxvol, clusterfilename.c_str(), argc, argv);
  } else if (methodname.compare("multipagerank") == 0) {
    return multicluster_pagerank_indiv(graphfilename,
      maxvol, clusterfilename.c_str(), argc, argv);
  } else {
    fprintf(stderr,"unknown method name: %s\n", methodname.c_str());
  }
}
