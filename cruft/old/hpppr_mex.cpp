/**
 * @file hpppr_mex.cpp
 * @author Kyle Kloster, David F. Gleich
 */

/**
 *
 * This file implements (I - alpha*P)^-1 *b via an approximate matrix-vector method
 * that uses a heap to determine the largest magnitude entries
 * of the iterative vector, then using only those nonzeros
 * to do the matvec.
 *
 */

#include "mex.h"
#include <vector>
#include <assert.h>
#include <math.h>
#include "taydeg.hpp"
#include "minheapnoL.hpp" // include our heap functions

/**
 * @param n - sparse matrix size
 * @param cp - sparse matrix column pointer (CSC)
 * @param ari - sparse matrix row index (CSC)
 * @param a - sparse matrix values (CSC)
 * @param set - an array of column indices, representing seed nodes
 * @param N - number of terms of Taylor polynomial
 * @param clustersize - target number of nonzeros to have in answer, y (1 <= clustersize <= n)
 *          (called 'maxnnz' here to avoid confusing 'clustersize' with the size of 'set')
 * @param y - the output vector (length n)
 */

void hpppr(const mwSize n, const mwIndex* cp, const mwIndex* ari, const double* a,
                const std::vector<mwIndex>& set, const int N, const double alpha,
                const mwIndex maxnnz, double* y)
{
    
    /**
     *
     *          INITILIAZATIONS
     *
     **/
    
    mwSize maxlistsize = n;
	double valind,minval;
    mwIndex hsize = 0;
    mwIndex ind, Tind, ri, arinzi, listind, nzi, bnzi;
    mwSize listsize = 0; ind = 0; nzi = 0; bnzi = 0; listind = 0;
    
    // allocate data
    
    std::vector<mwIndex> Tvec(maxnnz,0);
    std::vector<double> Vvec(maxnnz,0);
    mwIndex *T = &Tvec[0];
    double *v = &Vvec[0];
	
    std::vector<mwIndex> listvec(maxlistsize,0);
    mwIndex *list = &listvec[0];
    
    
    /*************************************************
     *
     *            HORNER'S RULE ITERATES
     *
     * Each Horner iterate consists of two steps:
     *
     * (1) HEAP into v
     *
     *  Filter the entries of y through a heap,
     *  sending only the largest to v. Stop as soon
     *  as every entry of vec has been checked,
     *  or we've added enough entries to make
     *  the 1-norm > tol, whichever occurs first.
     *
     * ---- I've temporarily taken out the 'tol'
     *      feature until I can debug
     * ----     - Kyle, 10/10/13
     *
     *
     * (2) MATVEC y = A*v
     *
     *  - A is compressed sparse column (CSC)
     *  - for each entry in the heap, add v[T[i]]*A[:,T[i]] to y
     *  - do this for each Horner's Rule iterate, and then we're done.
     *
     * (3) Add b to y
     *
     *************************************************/
    for (int k=-1; k <= N-1 ; k++){
        listind = 0; hsize = 0; valind = 0.; minval = 0.; ri = 0;
        
        /*********************************************************
         * (1)            HEAP y into v
         *********************************************************/
        
        while ( listind < listsize ){
            ind = list[listind];
            valind = y[ind];
            if (valind > 0) { // if valind is nonzero, it might go in heap
                if ( hsize >= maxnnz ){  // if heap is full, we might have to delete
                    minval = y[T[0]];
                    // if valind is big enough, replace the min entry with valind
                    if ( minval < valind ){
                        // replace T[0] with ind, then heap_down
                        ri = T[0];
                        T[0] = ind;
                        y[ri] = 0; // drop that entry from y
                        heap_down(0, hsize, T, y);
                        minval = y[T[0]];
                    }
                    else y[ind] = 0; // y[ind] so small it's excluded from heap
                }
                else { // if heap not full, place y[ind] at back of heap, then heap_up
                    T[hsize] = ind;
                    hsize++;
                    heap_up(hsize-1, hsize, T, y);
                }
            }
            listind ++;
        }
        
        // Copy the nonzeros in the heap of y into v
        // and zero out y, so it can be set equal to
        // the matvec y = A*v
        for (ind = 0; ind < hsize ; ind++){
            Tind = T[ind];
            v[ind] = y[Tind];
            y[Tind] = 0;
        }
        
        
        
        /*********************************************************
         * (2)          ACTUAL MATVEC
         *
         * At this point, y[] = 0 entirely, and v contains
         * just the entries from y that survived the heap process,
         * i.e. the maxnnz largest entries that were in y.
         *********************************************************/
        listsize = 0;
        for (ind = 0; ind < hsize; ind++) {
            Tind = T[ind];
            valind = v[ind]*alpha;
            for ( nzi=cp[Tind]; nzi < cp[Tind+1]; ++nzi) {
                arinzi = ari[nzi];
                if (y[arinzi] == 0){
                    list[listsize] = arinzi;
                    listsize++;
                }
                y[arinzi] += a[nzi]*valind;
            }
            v[ind] = 0;
        }
        
        /*********************************************************
         * (3)          ADD b to y
         *********************************************************/
        for (nzi=0; nzi < set.size(); ++nzi) {
            ind = set[nzi];
            if (y[ind] == 0){
                list[listsize] = ind;
                listsize++;
            }
            y[ind] += 1;
        }
        
    }//Taylor terms in Horner's rule are complete
    
    return;
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
// [prvector] = hpppr_mex(A,set,clustersize,alpha,tol)
void mexFunction(int nargout, mxArray *pargout[],       // these are your outputs
                 int nargin, const mxArray *pargin[])   // these are your arguments
{
    // inputs
    // A - sparse n-by-n
    // set - array of indices of seed nodes
    // clustersize - integer, target size of cluster
    // alpha - scalar 0 < alpha < 1
    // tol - scalar, 0<tol<1, accuracy of function expmv(x)
    
    const mxArray* A = pargin[0];
    const mxArray* set = pargin[1];
    mwSize n = mxGetM(A);
    
    mwIndex clustersize = ceil(pow(n,.5));
    double alpha = .9;
    double tol = pow(10,-6);
    
    if (nargin == 5){
        clustersize = (mwIndex)mxGetScalar(pargin[2]);
        alpha = mxGetScalar(pargin[3]);
        tol = mxGetScalar(pargin[4]);
    }
    else if (nargin == 4){
        clustersize = (mwIndex)mxGetScalar(pargin[2]);
        alpha = mxGetScalar(pargin[3]);
    }
    else if (nargin == 3){
        clustersize = (mwIndex)mxGetScalar(pargin[2]);
    }
    
    int N = ceil(log(tol*(1-alpha))/log(alpha))-1;

    std::vector< mwIndex > seeds;
    copy_array_to_index_vector( set, seeds );

    clustersize = clustersize*ceil(pow(log(n),2))*seeds.size();
    
    pargout[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    // decode the sparse matrix
    mwIndex* cp = mxGetJc(A);
    mwIndex* ri = mxGetIr(A);
    double* a = mxGetPr(A);
    double* y = mxGetPr(pargout[0]);
    
    mxAssert( nargin < 2, "incorrect number of inputs");
    if(tol<=0) {
        mexErrMsgIdAndTxt( "MATLAB:hpexpmv_mex:invalidInputs",
                          "tol must be positive.");
    }else if(nargin>5) {
        mexErrMsgIdAndTxt( "MATLAB:hpexpmv_mex:excessInputs",
                          "Too many inputs.");
    }
    mxAssert( tol > 0, "tol must be bigger than 0");
    mxAssert( 1 <= clustersize && clustersize <= n, "we must have 1 <= clustersize <= n");
    mxAssert( 0 < alpha && alpha < 1, "we must have 0 < alpha < 1");
    
    hpppr(n, cp, ri, a, // sparse matrix
               seeds, // seed vector
               N, alpha, clustersize,// parameters
               y); // output product vector
    return;
}
