/** 
 * @file gexpmq_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 * This file implements an approximate gauss-southwell method using a Queue
 * instead of a heap to approximate the largest element for the truncated 
 * taylor series approximation for a column of the matrix exponential
 */

#include "mex.h"
#include <queue>
#include <vector>
#include <assert.h>
#include <math.h>

/** A replacement for std::queue<int> using a circular buffer array */
class array_queue {
    public:
    size_t max_size;
    std::vector<int> array;
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

/**
 * @param n - sparse matrix size
 * @param cp - sparse matrix column pointer (CSC)
 * @param ari - sparse matrix row index (CSC)
 * @param a - sparse matrix values (CSC)
 * @param set - the set of seed nodes (entries of the preference vector that = 1)
 * @param alpha - the "teleportation constant" in pagerank (0 < alpha < 1)
 * @param tol - the stopping tolerance (0 < tol < Inf)
 * @param maxsteps - the maximum number of steps to take (1 <= maxsteps < Inf)
 * @param y - the output vector (length n)
 * @param nsteps - the number of output steps (length 1)
 */
void gsqppr(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a,
            const std::vector<mwIndex>& set, const double alpha,  const double tol, const mwIndex maxsteps,
            double* y, double *nsteps, double *npushes)
{
    //mexPrintf("Input n=%i N=%i c=%i tol=%i maxsteps=%i\n", n, alpha, c, tol, maxsteps);
    double sumresid = 0.;
    double toln = tol/n;
    // allocate data 
    std::vector<double> rvec(n,0.);
    double *r = &rvec[0];
    
    std::queue<mwIndex> Q;
    
    // set the initial residual, add to the heap, and update
    for (size_t i=0; i<set.size(); ++i) {
        r[set[i]]=1;
        Q.push(set[i]);
    }
    sumresid += (double)set.size();
    
    mwIndex npush = 0;
    *nsteps = (double)maxsteps; // set the default, which we change on early exit
    
    
    for (mwIndex iter = 0; iter < maxsteps; ++iter) {
        /* STEP 1: pop top element off of heap
         *  * get index i from T
         *  * add r(i) to y(i)
         *  * set r(i) to zero (update sumresid)
         * STEP 2: get i^th column from A
         *  * get neighbors of ith node
         *  * add as a column to r, scaled by alpha, and update heap
         *  *  (update sumresid)
         * Check for convegence!
        */
        
        // STEP 1: pop top element off of heap
        mwIndex ri = Q.front();
        Q.pop();
        
        double ri_val = r[ri];
        
        // update yi
        y[ri] += ri_val;
        
        // update r, no need to update heap here 
        r[ri] = 0;
        sumresid -= ri_val;
        double ri_valalpha = ri_val*alpha;
        
            // we add the column of A to the residual, scaled
            for (mwIndex nzi=cp[ri]; nzi < cp[ri+1]; ++nzi) {
                mwIndex v = ari[nzi];
                double ajv = a[nzi];
                double reold = r[v];
                r[v] += ajv*ri_valalpha;
                sumresid += ajv*ri_valalpha;
                if (r[v] > toln) {
                    if (reold < toln) {
                        Q.push(v);
                    }
                }
            }
            npush+=cp[ri+1]-cp[ri];
        //        if (sumresid < tol || Q.size() == 0 || sumsol > -tol ) {
        if (sumresid < tol || Q.size() == 0) {
            *nsteps = (double)iter;
            break;
        }
    } // end of for loop

    *npushes = (double)npush;
    return; // because we "break" out of for loop
}



void copy_array_to_index_vector(const mxArray* v, std::vector<mwIndex>& vec)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);
    
    vec.resize(n);
    
    for (size_t i=0; i<n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1 && n>=elem, "input 'set' contains entries violating n>= set(i) >= 1");
        vec[i] = (mwIndex)elem - 1;
    }
}



// USAGE
// [pagerankvector, nsteps, npushes] = gsqppr_mex(A,set,alpha,tol,maxsteps)
void mexFunction(
  int nargout, mxArray *pargout[],       // these are your outputs
  int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // set - vector of indices of starting nodes, 1 <= set(i) <= n
    // alpha - scalar 0 < alpha < 1
    // tol - double value, 0 < tol < Inf
    // maxsteps - integer scalar max-steps
    
    const mxArray* A = pargin[0];
    mwSize n = mxGetM(A);
    const mxArray* set = pargin[1];
    
    double alpha = 0.9;
    if (nargin >= 3){
        alpha = mxGetScalar(pargin[2]);
    }
    
    double tol = pow(10,-7);
    if (nargin >=4){
        tol = mxGetScalar(pargin[3]);
    }
//    tol = tol*(1-alpha); // need ||r|| < tol(1-alpha) for error < tol

    mwIndex maxsteps = 100*n;
    if (nargin == 5){
        maxsteps = (mwIndex)mxGetScalar(pargin[4]);
    }
    
    std::vector< mwIndex > cluster;
    copy_array_to_index_vector( set, cluster );

    pargout[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    pargout[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    pargout[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    // decode the sparse matrix
    mwIndex* cp = mxGetJc(A);
    mwIndex* ri = mxGetIr(A);
    double* a = mxGetPr(A);
    
    double* y = mxGetPr(pargout[0]);
    double* nsteps = mxGetPr(pargout[1]);
    double* npushes = mxGetPr(pargout[2]);
    
    // note: these 'Assert' error messages won't work
    // unless the mex code is compiled with the '-g' option
    
    mxAssert(alpha > 0 && alpha < 1, "alpha must satisfy 0 < alpha < 1");
    mxAssert(tol > 0 && tol <= 1, "tol must be 0 < tol <= 1");
    mxAssert(maxsteps > 0, "we must have maxsteps > 0");
    
    gsqppr(n, cp, ri, a, // sparse matrix
           cluster, alpha, tol, maxsteps, // parameters
           y, nsteps, npushes);
    
}
    
    
  