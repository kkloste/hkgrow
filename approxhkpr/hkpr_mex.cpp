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
#include <random>
#include <ctime>
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

int debugflag = 0;


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

/*****
 *
 *          above:  DATA STRUCTURES
 *
 *
 *
 *          below:  CLUSTERING FUNCTIONS
 *
 ****/

long unsigned int poisson(long unsigned int lambda){
}

mwIndex random_walk(sparserow * G, mwIndex K,  mwIndex cur_node){
    mwIndex next_node, nzi, v;
        std::default_random_engine generator;
    for(mwIndex steps = 1; steps <= K; steps++){
        std::uniform_int_distribution<int> distribution(0,sr_degree(G,cur_node) - 1);
        next_node = distribution(generator);
        nzi=G->ai[cur_node];
        cur_node = G->aj[nzi + next_node];
    }

    return cur_node;
}

int approx_hkpr(sparserow * G, sparsevec& y, const double alphat, std::vector<mwIndex> seedvec, const double tol, int debugflag){
    long unsigned int lambda = (long unsigned int)ceil(alphat);
    mwIndex n = G->n;
    long unsigned int v,k;
    double r = (16.0/pow(tol,2))*log(n);
    long unsigned int K = (long unsigned int)ceil(log(1.0/tol)/log(log(1.0/tol)));
    std::mt19937 mrand(std::time(0));  // seed however you want
    std::poisson_distribution<long unsigned int> d(lambda);
    size_t seedsize = seedvec.size(); // first pick seed node uniformly at random from seedvec.
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, seedsize - 1);
    
    mwIndex cur_node;

if (debugflag>=1){ mexPrintf("r=%f K=%i \n", r,K); mexEvalString("drawnow"); }

    for (mwIndex iter = 1; iter <= (int)ceil(r); iter++){
        k = d(mrand); // determines length of walk
if (debugflag>=2){ mexPrintf("k=%i ", k); mexEvalString("drawnow"); }
        k = std::min(k,K);
        cur_node = seedvec[distribution(generator)];
        v = random_walk(G, k, cur_node);
if (debugflag>=2){ mexPrintf("v=%i ", v); mexEvalString("drawnow"); }
        y.map[v] += 1;
    }

    return r;
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
    std::vector<int> volume(prpairs.size());
    std::vector<int> cutsize(prpairs.size());
    
    size_t i=0;
    tr1ns::unordered_map<int,size_t> rank;
    for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
         it!=itend; ++it, ++i) {
        rank[it->first] = i;
    }
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
 * @param alpha the value of alpha in the heatkernelPageRank computation
 * @param target_vol the approximate number of edges in the cluster
 * @param p the heatkernelpagerank vector
 * @param r the residual vector
 * @param a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */
template <class Queue>
int hypercluster_heatkernel_multiple(sparserow* G,
                                     const std::vector<mwIndex>& set, double alpha, double target_vol,
                                     sparsevec& p, sparsevec &r, Queue& q,
                                     std::vector<mwIndex>& cluster, local_hkpr_stats *stats, double tol)
{
    // reset data
    p.map.clear();
    r.map.clear();
    q.empty();
    
    assert(target_vol > 0);
    //    assert(alpha < 1.0); assert(alpha > 0.0);
    // this is commented out because alpha will be >=1 for exp(tG), in contrast with pagerank
    
    
    //r.map[start] = 1.0;
    size_t maxdeg = 0;
    for (size_t i=0; i<set.size(); ++i) { //populate r with indices of "set"
        assert(set[i] >= 0); assert(set[i] < G->n); // assert that "set" contains indices i: 1<=i<=n
        size_t setideg = sr_degree(G,set[i]);
        //        r.map[set[i]] = (1./(double)(set.size()))/(double)setideg;
        r.map[set[i]] = 1.;
        maxdeg = std::max(maxdeg, setideg);
    }
    
    
    /*
     //double pr_eps = 1.0/std::max((double)sr_degree(G,start)*(double)target_vol, 100.0);
     //double pr_eps = std::min(1.0/std::max(10.*target_vol, 100.0),
     //1./(double)(set.size()*maxdeg + 1));
     double pr_eps = 1.0/std::max(10.0*target_vol, 100.0);
     if (stats) { stats->eps = pr_eps; }
     
     // calculate an integer number of max_push_count
     double max_push_count = 1./(pr_eps*(1.-alpha));
     max_push_count = std::min(max_push_count, 0.5*(double)std::numeric_limits<int>::max());
     */
    
    //      instead of pr_eps, use tol = 1e-8
    //      instead of '(int)max_pust_count', use ceil(pow(G->n,1.5))
    
    if (debugflag >= 1){ mexPrintf("approx_hkpr: alphat=%f tol=%f \n", alpha, tol); mexEvalString("drawnow");}

    int nsteps = approx_hkpr(G, p, alpha, set, tol, debugflag);
    
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




/**
 *          HKGROW
 *  G = sparserow version of input matrix A
 *  set = vector of seed nodes
 **/
void hkgrow(sparserow* G, std::vector<mwIndex>& seeds, double alpha,
            double targetvol, double* fcond, double* fcut,
            double* fvol, double tol)
{
    sparsevec p, r;
    std::queue<mwIndex> q;
    local_hkpr_stats stats;
    std::vector<mwIndex> bestclus;
    if (debugflag >= 1){ mexPrintf("hkpr_mex: call to hypercluster_heatkernel() start \n"); mexEvalString("drawnow");}
    hypercluster_heatkernel_multiple(G, seeds, alpha, targetvol,
                                     p, r, q, bestclus, &stats, tol);
    if (debugflag >= 1){ mexPrintf("hkpr_mex: call to hypercluster_heatkernel() DONE \n"); mexEvalString("drawnow");}
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
// [bestset,cond,cut,vol] = hkpr_mex(A,set,targetvol,alpha,tol,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    debugflag = (int)mxGetScalar(prhs[5]);
    if (debugflag >= 1){ mexPrintf("hkpr_mex: preprocessing start: \n");mexEvalString("drawnow");}
    
    
    mxAssert(nrhs > 2 && nrhs < 7, "2-6 inputs required.");
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    mxAssert(mxIsSparse(mat), "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");
    
    mxArray* cond = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(1,1,mxREAL);
    
    if (debugflag >= 1){mexPrintf("hkpr_mex: declared some input/outputs: \n");    mexEvalString("drawnow");}
    
    if (nlhs > 1) { plhs[1] = cond; }
    if (nlhs > 2) { plhs[2] = cut; }
    if (nlhs > 3) { plhs[3] = vol; }
    
    mxAssert(nlhs <= 4, "Too many output arguments");
    
    double tol = pow(10,-5);
    double alpha = 1.;
    
    
    if (debugflag >= 1){mexPrintf("hkpr_mex: declared more input/outputs: \n");    mexEvalString("drawnow");}
    
    if (nrhs >= 4) {
        alpha = mxGetScalar(prhs[3]);
        tol = mxGetScalar(prhs[4]);
    }
    
    // use a strange sentinal
    double targetvol = 1000.;
    if (nrhs >= 3) {
        targetvol = mxGetScalar(prhs[2]);
    }
    
    
    if (debugflag >= 1){mexPrintf("hkpr_mex: input/outputs 3 : \n");    mexEvalString("drawnow");}
    
    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);
    
    if (debugflag >= 1){mexPrintf("hkpr_mex: input/outputs 4 : \n");    mexEvalString("drawnow");}
    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    
    if (debugflag >= 1){mexPrintf("hkpr_mex: preprocessing end: \n"); mexEvalString("drawnow");}
    
    hkgrow(&r, seeds, alpha, targetvol,
           mxGetPr(cond), mxGetPr(cut), mxGetPr(vol), tol);
    if (debugflag >= 1){mexPrintf("hkpr_mex: call to hkgrow() done \n"); mexEvalString("drawnow");}
    
    if (nlhs > 0) {
        mxArray* cassign = mxCreateDoubleMatrix(seeds.size(),1,mxREAL);
        plhs[0] = cassign;
        
        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<seeds.size(); ++i) {
            ci[i] = (double)(seeds[i] + 1);
        }
    }
}