/**
 * @file greedyclustergrow.cc
 * Greedily grow a cluster to minimize conductance
 */

#include <vector>
#include <queue>
#include <utility> // for pair sorting

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


double sparse_array_volume(const mxArray* mat)
{
    mwIndex *aj = mxGetIr(mat), *ai = mxGetJc(mat);
    mwSize n = mxGetM(mat);
    double vol = 0;
    mwIndex last = aj[n+1];
    for (mwIndex i=0; i<n; ++i) {
        for (mwIndex nzi=ai[i]; nzi<ai[i+1]; ++nzi) {
            if (i != aj[i]) {
                vol += 1.;
            }
        }
    }

    return vol;
}

void greedy_cluster_grow(const mxArray* mat, double* fcond, double* fcut,
    double* fvol, std::vector<size_t >& cluster, 
    double Gfullvol, double deg_pseudo_count)
{
    // we only handle symmetric matrices
    mwIndex *aj = mxGetIr(mat), *ai = mxGetJc(mat);
    mwSize n = mxGetM(mat);

    tr1ns::unordered_set<int> ind;
    tr1ns::unordered_map<mwIndex, mwIndex> bedges;  
                                        // bound the number of boundary
                                        // edges incident on a vertex

    // index the starting set
    double avgdeg = 0.;
    for (size_t i=0; i<cluster.size(); ++i) {
        mwIndex v = cluster[i];
        ind.insert(v);
        avgdeg += ai[v+1]-ai[v]; // loops don't matter for this estimate
    }
    avgdeg /= (double)cluster.size();

    // use weird sentinal values
    if (deg_pseudo_count >= -100.*n && deg_pseudo_count <= 100.*n) {
        avgdeg = deg_pseudo_count;
    }

    // compute the initial cut/volume
    // along with the boundary queue
    double cut = 0;
    double vol = 0;
    std::priority_queue< std::pair<double,mwIndex> > queue;

    for (size_t i=0; i<cluster.size(); ++i) {
        mwIndex v = cluster[i];
        for (mwIndex nzi=ai[v]; nzi<ai[v+1]; ++nzi) {   
            mwIndex x = aj[nzi];
            if (v == x) { 
                continue;
            }
            vol += 1.;

            if (ind.count(x) == 0) {
                cut += 1.;
                bedges[x] += 1;
                // implement a really silly updated priority queue
                queue.push(std::make_pair( 
                    (double)bedges[x]/(avgdeg+(double)ai[x+1]-(double)ai[x]), x ));
            } 
        }
    }

    // now greedily grow the cluster
    size_t nextra = 0; // the number of nodes in the queue we need to zero
    tr1ns::unordered_set<mwIndex> queue_used;

    while (queue.size() > 0) {
        std::pair< double, mwIndex > head = queue.top();
        queue.pop();
        if (queue_used.count(head.second) > 0) { continue; }

        mwIndex w = head.second;

        // check if this is the most recent score or not
        if (head.first != (double)bedges[w]/(avgdeg+(double)ai[w+1]-(double)ai[w])) {
            //mexPrintf("%i: skipped because %f != %f\n", v, head.first, 
                //(double)bedges[w]/(myd+(double)ai[w+1]-(double)ai[w]));
            continue; 
        }

        queue_used.insert(w);

        bedges[w] = 0.;

        double newcut = cut;
        double newvol = vol;

        // see how cut and vol will change if we add w to the cluster
        for (mwIndex nzi2=ai[w]; nzi2<ai[w+1]; ++nzi2) {
            mwIndex x = aj[nzi2];
            if (x == w) { continue; }
            newvol += 1.;
            if (ind.count(x) > 0) {
                newcut -= 1.;
            } else {
                newcut += 1.;
                bedges[x] += 1;
                // implement a really silly updated priority queue
                queue.push(std::make_pair( 
                    (double)bedges[x]/(avgdeg+(double)ai[x+1]-(double)ai[x]), x ));
            }
        }

        double newcond = newcut/std::min(newvol,Gfullvol-newvol);
        double prevcond = cut/std::min(vol,Gfullvol-vol);

        if (newcond < prevcond) {
            //mexPrintf("%i: added %i dcond=%f\n", v, w, newcond-prevcond);
            ind.insert(w);
            cluster.push_back(w);
        } else {
            //mexPrintf("%i: ended at %i dcond=%f\n", v, w, newcond-prevcond);
            break;
        }

        cut = newcut;
        vol = newvol;
    }

    *fcut = cut;
    *fvol = vol;
    *fcond = cut/std::min(vol,Gfullvol-vol);
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

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    mxAssert(nrhs > 2 && nrhs < 5, "2-4 inputs required.");

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

    std::vector< size_t > cluster;
    copy_array_to_index_vector( set, cluster );

    double Gvol = 0;
    if (nrhs >= 3) {
        Gvol = mxGetScalar(prhs[2]);
    } else {
        Gvol = sparse_array_volume(mat);
    }

    // use a strange sentinal
    double deg_pseudo = -1.*mxGetM(mat)*100. - 100;
    if (nrhs >= 4) {
        deg_pseudo = mxGetScalar(prhs[3]);
    }

    greedy_cluster_grow(mat, 
        mxGetPr(cond), mxGetPr(cut), mxGetPr(vol),
        cluster, Gvol, deg_pseudo);

    if (nlhs > 0) { 
        mxArray* cassign = mxCreateDoubleMatrix(cluster.size(),1,mxREAL);
        plhs[0] = cassign;

        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<cluster.size(); ++i) {
            ci[i] = (double)(cluster[i] + 1);
        }
    }
}