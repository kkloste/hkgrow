/**
 * @file hkclus_mex.cc
 * Implement a personal heat kernel pagerank clustering scheme.
 *
 * mex hkclus_mex.cc CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims
 *
 *
 *
 *
 */


#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include <math.h>

#ifdef __APPLE__
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#define tr1ns std::tr1
#else
#include <unordered_set>
#include <unordered_map>
#define __STDC_UTF_16__ 1
#define tr1ns std
#endif

#include <mex.h>

#define DEBUGPRINT(x) do { if (debugflag) { \
                            mexPrintf x; mexEvalString("drawnow"); } \
                      } while (0)

int debugflag = 0;

struct sparsevec {
    typedef tr1ns::unordered_map<mwIndex,double> map_type;
    map_type map;
    /** Get an element and provide a default value when it doesn't exist
     * This command does not insert the element into the vector
     */
    double get(mwIndex index, double default_value=0.0) {
        map_type::iterator it = map.find(index);
        if (it == map.end()) {
            return default_value;
        } else {
            return it->second;
        }
    }
    
    /** Compute the sum of all the elements
     * Implements compensated summation
     */
    double sum() {
        double s=0.;
        for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
            s += it->second;
        }
        return s;
    }
    
    /** Compute the max of the element values
     * This operation returns the first element if the vector is empty.
     */
    mwIndex max_index() {
        mwIndex index=0;
        double maxval=std::numeric_limits<double>::min();
        for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
            if (it->second>maxval) { maxval = it->second; index = it->first; }
        }
        return index;
    }
};

struct sparserow {
    mwSize n, m;
    mwIndex *ai;
    mwIndex *aj;
    double *a;
};


mwIndex sr_degree(sparserow *s, mwIndex u) {
    return (s->ai[u+1] - s->ai[u]);
}


int taylordegree(const double t, const double eps) {
    double eps_exp_t = eps*exp(t);
    double error = exp(t)-1;
    double last = 1.;
    double k = 0.;
    while(error > eps_exp_t){
        k = k + 1.;
        last = (last*t)/k;
        error = error - last;
    }
    return std::max((int)k, (int)3);
}

/*****
 *
 *          above:  DATA STRUCTURES
 *
 *
 *
 *          below:  CLUSTERING FUNCTIONS
 *
 ****/

/**
 *
 *  gsqexpmseed inputs:
 *      G   -   adjacency matrix of an undirected graph
 *      set -   seed vector: the indices of a seed set of vertices
 *              around which cluster forms.
 *              Rather than normalize 'set' (by setting
 *                  set[i] = 1/set.size(); )
 *              we instead multiply eps by set.size().
 *  output: 
 *      y = exp(-t(I-P)) * set
 *              with infinity-norm accuracy of eps
 *              in the degree weighted norm
 *  parameters:
 *      t   - the value of t
 *      eps - the accuracy
 *      max_push_count - the total number of steps to run
 *      Q - the queue data structure
 */
template <class Queue>
int gsqexpmseed(sparserow * G, sparsevec& set, sparsevec& y,
                const double t, const double eps,
                const mwIndex max_push_count, Queue& Q)
{
    DEBUGPRINT(("gsqexpmseed interior: t=%f eps=%f \n", t, eps)); 
    mwIndex n = G->n;
    mwIndex N = (mwIndex)taylordegree(t, eps);
    DEBUGPRINT(("gsqexpmedseed: n=%i N=%i \n", n, N));
    
    // initialize the weights for the different residual partitions
    // r(i,j) > d(i)*psi_1(t)*eps / (N*psi_j(t))
    //  since each coefficient but d(i) stays the same,
    //  we combine all coefficients except d(i)
    //  into the vector "pushcoeff"
    std::vector<double> psivec(N+1,0.);
    psivec[N] = 1;
    for (int k = 1; k < N ; k++){
        psivec[N-k] = psivec[N-k+1]*t/(double)(N-k+1) + 1;
    } // psivec[k] = psi_k(t)
    double eps_exp_t = eps*psivec[1];
    std::vector<double> pushcoeff(N+1,0.);
    pushcoeff[1]=eps/(double)N;
    for (int k = 2; k <= N ; k++){
        pushcoeff[k] = pushcoeff[k-1]*(psivec[k-1]/psivec[k]);
    }
    pushcoeff[0]=0;
    
    mwIndex ri = 0;
    mwIndex npush = 0;
    double rij = 0;
    // allocate data
    sparsevec rvec;

    // i is the node index, j is the "step"
    #define rentry(i,j) ((i)+(j)*n)
    
    // set the initial residual, add to the queue
    for (sparsevec::map_type::iterator it=set.map.begin(),itend=set.map.end(); 
         it!=itend;++it) {
        ri = it->first;
        it->second = 1./(double)sr_degree(G,ri);
    }
    double scalefactor = set.sum();
    for (sparsevec::map_type::iterator it=set.map.begin(),itend=set.map.end();
         it!=itend;++it) {
        ri = it->first;
        rij = scalefactor*(double)it->second;
        rvec.map[rentry(ri,0)]+=rij;
        Q.push(rentry(ri,0));
    }
    
    
    while (npush < max_push_count) {
        // STEP 1: pop top element off of heap
        ri = Q.front();
        Q.pop();
        // decode incides i,j
        mwIndex i = ri%n;
        mwIndex j = ri/n;
        
        double degofi = (double)sr_degree(G,i);
        double kappai = degofi*pushcoeff[j];
        rij = rvec.map[ri];
//        
        // update yi
        y.map[i] += rij;
        
        // update r, no need to update heap here
        rvec.map[ri] = 0; 
        
        double rijs = t*rij/(double)(j+1);
        double ajv = 1./degofi;
        double update = rijs*ajv;
        
        if (j == N-1) {
            // this is the terminal case, and so we add the column of A
            // directly to the solution vector y
            for (mwIndex nzi=G->ai[i]; nzi < G->ai[i+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                y.map[v] += ajv*rijs;
            }
            npush += degofi;
        }
        else {
            // this is the interior case, and so we add the column of A
            // to the residual at the next time step.
            for (mwIndex nzi=G->ai[i]; nzi < G->ai[i+1]; ++nzi) {
                mwIndex v = G->aj[nzi];
                mwIndex re = rentry(v,j+1);
                double reold = rvec.get(re);
                double renew = reold + update;
                double dv = sr_degree(G,v);
                rvec.map[re] = renew;
                if (renew >= dv*pushcoeff[j+1] && reold < dv*pushcoeff[j+1]) {
                    Q.push(re);
                }
            }
            npush+=degofi;
        }
        // terminate when Q is empty, i.e. we've pushed all r(i,j) > eps*psi_1(t)*d(i)/(N*psi_j(t))
        if ( Q.size() == 0) { return npush; }
    }//end 'while'
    return (npush);
}


struct greater2nd {
    template <typename P> bool operator() (const P& p1, const P& p2) {
        return p1.second > p2.second;
    }
};

void cluster_from_sweep(sparserow* G, sparsevec& p, 
      std::vector<mwIndex>& cluster, double *outcond, double* outvolume,
      double *outcut)
{
  // now we have to do the sweep over p in sorted order by value
  typedef std::vector< std::pair<int, double> > vertex_prob_type;
  vertex_prob_type prpairs(p.map.begin(), p.map.end());
  std::sort(prpairs.begin(), prpairs.end(), greater2nd());

  // compute cutsize, volume, and conductance
  std::vector<double> conductance(prpairs.size());
  std::vector<mwIndex> volume(prpairs.size());
  std::vector<mwIndex> cutsize(prpairs.size());

  size_t i=0;
  tr1ns::unordered_map<int,size_t> rank;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    rank[it->first] = i;
  }
  //printf("support=%i\n",prpairs.size());
  mwIndex total_degree = G->ai[G->m];
  mwIndex curcutsize = 0;
  mwIndex curvolume = 0;
  i=0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    mwIndex v = it->first;
    mwIndex deg = G->ai[v+1]-G->ai[v];
    mwIndex change = deg;
    for (mwIndex nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
      mwIndex nbr = G->aj[nzi];
      if (rank.count(nbr) > 0) {
        if (rank[nbr] < rank[v]) {
          change -= 2;
        }
      }
    }
    curcutsize += change;
    //if (curvolume + deg > target_vol) {
      //break;
    //}
    curvolume += deg;
    volume[i] = curvolume;
    cutsize[i] = curcutsize;
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
  if (outcut) { *outcut = cutsize[mincondind]; }
}

struct local_hkpr_stats {
    double conductance;
    double volume;
    double support;
    double steps;
    double eps;
    double cut;
};

/** Cluster will contain a list of all the vertices in the cluster
 * @param set the set of starting vertices to use
 * @param t the value of t in the heatkernelPageRank computation
 * @param eps the solution tolerance eps
 * @param p the heatkernelpagerank vector
 * @param r the residual vector
 * @param a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */
template <class Queue>
int hypercluster_heatkernel_multiple(sparserow* G,
        const std::vector<mwIndex>& set, double t, double eps,
        sparsevec& p, sparsevec &r, Queue& q,
        std::vector<mwIndex>& cluster, local_hkpr_stats *stats)
{
    // reset data
    p.map.clear();
    r.map.clear();
    q.empty();
    DEBUGPRINT(("beginning of hypercluster \n"));

    size_t maxdeg = 0;
    for (size_t i=0; i<set.size(); ++i) { //populate r with indices of "set"
        assert(set[i] >= 0); assert(set[i] < G->n); // assert that "set" contains indices i: 1<=i<=n
        size_t setideg = sr_degree(G,set[i]);
        r.map[set[i]] = 1./(double)(set.size()); // r is normalized to be stochastic
//    DEBUGPRINT(("i = %i \t set[i] = %i \t setideg = %i \n", i, set[i], setideg));
        maxdeg = std::max(maxdeg, setideg);
    }
    
    DEBUGPRINT(("at last, gsqexpm: t=%f eps=%f \n", t, eps));
    
    int nsteps = gsqexpmseed(G, r, p, t, eps, ceil(pow(G->n,1.5)), q);
/**
 *      **********
 *
 *        ***       GSQEXPMSEED       is called         ***
 *
 *      **********
 */

    if (nsteps == 0) {
        p = r; // just copy over the residual
    }
    int support = r.map.size();
    if (stats) { stats->steps = nsteps; }
    if (stats) { stats->support = support; }
    
 // scale the probablities by their degree
    for (sparsevec::map_type::iterator it=p.map.begin(),itend=p.map.end();
         it!=itend;++it) {
        it->second *= (1.0/(double)std::max(sr_degree(G,it->first),(mwIndex)1));
    }
    
    double *outcond = NULL;
    double *outvolume = NULL;
    double *outcut = NULL;
    if (stats) { outcond = &stats->conductance; }
    if (stats) { outvolume = &stats->volume; }
    if (stats) { outcut = &stats->cut; }
    cluster_from_sweep(G, p, cluster, outcond, outvolume, outcut);
    return (0);
}

/** Grow a set of seeds via the heat-kernel.
 *
 * @param G sparserow version of input matrix A
 * @param seeds a vector of input seeds seeds (index 0, N-1), and then
 *          updated to have the final solution nodes as well.
 * @param t the value of t in the heat-kernel
 * @param eps the solution tolerance epsilon
 * @param fcond the final conductance score of the set.
 * @param fcut the final cut score of the set
 * @param fvol the final volume score of the set
 */
void hkgrow(sparserow* G, std::vector<mwIndex>& seeds, double t,
             double eps, double* fcond, double* fcut,
             double* fvol)
{
    sparsevec p, r;
    std::queue<mwIndex> q;
    local_hkpr_stats stats;
    std::vector<mwIndex> bestclus;
    DEBUGPRINT(("hkgrow_mex: call to hypercluster_heatkernel() start\n"));
    hypercluster_heatkernel_multiple(G, seeds, t, eps,
                                   p, r, q, bestclus, &stats);
    DEBUGPRINT(("hkgrow_mex: call to hypercluster_heatkernel() DONE\n"));
    seeds = bestclus;
    *fcond = stats.conductance;
    *fcut = stats.cut;
    *fvol = stats.volume;
}

void copy_array_to_index_vector(const mxArray* v, std::vector<mwIndex>& vec)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);
    
    vec.resize(n);
    
    for (size_t i=0; i<n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1, "Only positive integer elements allowed");
        vec[i] = (mwIndex)elem - 1;
    }
}


// USAGE
// [bestset,cond,cut,vol] = hkgrow_mex(A,set,t,eps,debugflag)
// Note that targetvol is currently ignored
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 5) { 
        mexErrMsgIdAndTxt("hkgrow_mex:notEnoughArguments", 
            "hkgrow_mex needs six arguments not %i", nrhs);
    }
    debugflag = (int)mxGetScalar(prhs[4]);
    DEBUGPRINT(("hkgrow_mex: preprocessing start: \n"));
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    mxAssert(mxIsSparse(mat), "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");
    
    mxArray* cond = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(1,1,mxREAL);
    
    if (nlhs > 1) { plhs[1] = cond; }
    if (nlhs > 2) { plhs[2] = cut; }
    if (nlhs > 3) { plhs[3] = vol; }
    
    mxAssert(nlhs <= 4, "Too many output arguments");
    
    double eps = pow(10,-3);
    double t = 15.;
    
    if (nrhs >= 4) {
        t = mxGetScalar(prhs[2]);
        eps = mxGetScalar(prhs[3]);
    }
    
    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);
    
    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );

    DEBUGPRINT(("hkgrow_mex: preprocessing end: \n"));

    hkgrow(&r, seeds, t, eps, 
            mxGetPr(cond), mxGetPr(cut), mxGetPr(vol) );
    
    DEBUGPRINT(("hkgrow_mex: call to hkgrow() done\n"));
    
    if (nlhs > 0) { 
        mxArray* cassign = mxCreateDoubleMatrix(seeds.size(),1,mxREAL);
        plhs[0] = cassign;
        
        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<seeds.size(); ++i) {
            ci[i] = (double)(seeds[i] + 1);
        }
    }
}