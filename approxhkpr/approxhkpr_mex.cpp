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
 *          below:  APPROX_HKPR
 *
 ****/

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

int approx_hkpr(sparserow * G, double* y, const double alphat, std::vector<mwIndex> seedvec, const double tol, int debugflag){
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
        y[v] += 1;
    }
    
    return r;
}


struct greater2nd {
    template <typename P> bool operator() (const P& p1, const P& p2) {
        return p1.second > p2.second;
    }
};



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
// [x] = approxhkpr_mex(A,set,alpha,tol,debugflag)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    debugflag = (int)mxGetScalar(prhs[4]);
    if (debugflag >= 1){ mexPrintf("hkpr_mex: preprocessing start: \n");mexEvalString("drawnow");}
    
    
    mxAssert(nrhs > 1 && nrhs < 6, "2-5 inputs required.");
    
    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];
    
    mxAssert(mxIsSparse(mat), "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");
    
    mxAssert(nlhs <= 2, "Too many output arguments");
    mxAssert(nlhs > 0, "Not enough output arguments");
    double tol = pow(10,-1.5);
    double alpha = 1.;
    
    if (nrhs >= 4) { tol = mxGetScalar(prhs[4]); }
    if (nrhs >= 3) { alpha = mxGetScalar(prhs[3]); }
    
    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);
    
    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );
    
    if (debugflag >= 1){mexPrintf("hkpr_mex: preprocessing end: \n"); mexEvalString("drawnow");}

    plhs[0] = mxCreateDoubleMatrix(r.n,1,mxREAL);
    double* y = mxGetPr(plhs[0]);
    int npushes = approx_hkpr(&r, y, alpha, seeds, tol, debugflag);
    if (nlhs == 2) {
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        double* numpushes = mxGetPr(plhs[1]);
        numpushes[0] = (double)npushes;
    }
}