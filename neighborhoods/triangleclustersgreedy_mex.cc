/**
 * @file triangleclustersgreedy_mex.cc
 * Implement clustering coefficient calculations in C++ with
 * a greedy scheme for growing the clusters.
 */

#include <mex.h>

#include <vector>
#include <queue>
#include <utility> // for pair sorting

void triangle_clusters_greedy(const mxArray* mat, mxArray* cond, mxArray* cut,
    mxArray* vol, mxArray *s, mxArray* cc, mxArray* t, 
    std::vector<std::pair<size_t,size_t> >& clusters, bool verbose)
{
    // we only handle symmetric matrices
    mwIndex *aj = mxGetIr(mat), *ai = mxGetJc(mat);
    mwSize n = mxGetM(mat);

    double Gvol = 0.;

    std::vector<int> ind(n, 0);
    std::vector<double> priority(n,0.0);
    std::vector<mwIndex> bedges(n,0);  // bound the number of boundary
                                       // edges incident on a vertex

    std::vector<int> queue_used(n,0);
    std::vector<mwIndex> queue_extra(n,0);

    double Gfullvol = (double)ai[n];

    for (mwIndex v = 0; v < n; ++v) {
        //mexPrintf("%i\n", v);
        if (verbose && v%100==0) {
            mexPrintf("%5.1f%% %7i/%7i\n",
                100.*(double)(v+1.)/(double)n, v, n);
        }
        // index the neighbors
        for (mwIndex nzi=ai[v]; nzi<ai[v+1]; ++nzi) {
            ind[aj[nzi]] = 1;
            clusters.push_back(std::make_pair(v, aj[nzi]));
        }
        ind[v] = 1;
        clusters.push_back(std::make_pair(v, v));

        double d = (double)(ai[v+1]-ai[v]);
        double myd = d; // just a version that doesn't update
        double curvol = 0.;
        double curcut = 0.;
        double curt = 0.;

        std::priority_queue< std::pair<double,mwIndex> > queue;

        // do two BFS steps
        for (mwIndex nzi=ai[v]; nzi<ai[v+1]; ++nzi) {
            mwIndex w = aj[nzi];
            if (v == w) { d-=1.; continue; }
            curvol += 1.;
            Gvol += 1.;
            for (mwIndex nzi2=ai[w]; nzi2<ai[w+1]; ++nzi2) {
                mwIndex x = aj[nzi2];
                if (x == w) { continue; }
                curvol += 1.;
                if (x==v) { continue; }
                if (ind[x]) {
                    curt += 1.;
                } else {
                    curcut += 1.;
                    bedges[x] += 1;
                    // implement a really silly updated priority queue
                    queue.push(std::make_pair( 
                        (double)bedges[x]/(myd+(double)ai[x+1]-(double)ai[x]), x ));
                }
            }
        }

        double csize = (double)d + 1.;
        size_t nextra = 0; // the number of nodes in the queue we need to zero

        while (queue.size() > 0) {
            std::pair< double, mwIndex > head = queue.top();
            queue.pop();
            if (queue_used[head.second]) { continue; }

            //mexPrintf("%i: popped unused %i, %f\n", v, head.second, head.first);

            mwIndex w = head.second;
            // check if this is the most recent score or not
            if (head.first != (double)bedges[w]/(myd+(double)ai[w+1]-(double)ai[w])) {
                //mexPrintf("%i: skipped because %f != %f\n", v, head.first, 
                    //(double)bedges[w]/(myd+(double)ai[w+1]-(double)ai[w]));
                continue; 
            }
            
            queue_used[w] = 1;
            queue_extra[nextra] = w;
            nextra += 1;
            bedges[w] = 0;

            double newcut = curcut;
            double newvol = curvol;
            
            // see how cut and vol will change if we add w to the cluster
            for (mwIndex nzi2=ai[w]; nzi2<ai[w+1]; ++nzi2) {
                mwIndex x = aj[nzi2];
                if (x == w) { continue; }
                newvol += 1.;
                if (ind[x]) {
                    newcut -= 1.;
                } else {
                    newcut += 1.;
                    bedges[x] += 1;
                    // implement a really silly updated priority queue
                    queue.push(std::make_pair( 
                        (double)bedges[x]/(myd+(double)ai[x+1]-(double)ai[x]), x ));
                }
            }
        
            double newcond = newcut/std::min(newvol,Gfullvol-newvol);
            double prevcond = curcut/std::min(curvol,Gfullvol-curvol);

            if (newcond < prevcond) {
                //mexPrintf("%i: added %i dcond=%f\n", v, w, newcond-prevcond);
                ind[w] = 1;
                clusters.push_back(std::make_pair(v, w));
                csize += 1.;
            } else {
                //mexPrintf("%i: ended at %i dcond=%f\n", v, w, newcond-prevcond);
                break;
            }

            curcut = newcut;
            curvol = newvol;
        }

        // finish popping vertices
        while (queue.size() > 0) {
            std::pair< double, mwIndex > head = queue.top();
            queue.pop();
            mwIndex w = head.second;
            bedges[w] = 0;
        }

        // assign the output
        mxGetPr(cut)[v] = curcut;
        mxGetPr(vol)[v] = curvol;
        if (d > 1.) {
            mxGetPr(cc)[v] = curt/(d*(d-1.));
        } else {
            mxGetPr(cc)[v] = 0.;
        }
        mxGetPr(s)[v] = std::min(n-csize,csize);
        mxGetPr(t)[v] = curt/2.;

        // clear the index
        for (mwIndex nzi=ai[v]; nzi<ai[v+1]; ++nzi) {
            ind[aj[nzi]] = 0;
        }
        for (size_t ei=0; ei<nextra; ++ei) {
            mwIndex w = queue_extra[ei];
            queue_used[w] = 0; // reset the used flag
            ind[w] = 0;
        }
        ind[v] = 0;

        for (size_t i=0; i<n; ++i) {
            mxAssert(queue_used[i] == 0, "queue_used not reset");
            mxAssert(ind[i] == 0, "ind not reset");
            mxAssert(bedges[i] == 0, "ind not reset");
        }
    }

    for (mwIndex v = 0; v < n; ++v) {
        mxGetPr(cond)[v] = mxGetPr(cut)[v]/std::min(mxGetPr(vol)[v],Gvol-mxGetPr(vol)[v]);
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    mxAssert(nrhs == 1 || nrhs==2, "One or two input required.");

    const mxArray* mat = prhs[0];

    mxAssert(mat, "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");

    mwSize n=mxGetM(mat);

    bool verbose = false; 
    if (nrhs == 2) {
        verbose = (bool)mxGetScalar(prhs[1]);
    }

    mxArray* cond = mxCreateDoubleMatrix(n,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(n,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(n,1,mxREAL);
    mxArray* s = mxCreateDoubleMatrix(n,1,mxREAL);
    mxArray* cc = mxCreateDoubleMatrix(n,1,mxREAL);
    mxArray* t = mxCreateDoubleMatrix(n,1,mxREAL);
    
    if (nlhs > 0) { plhs[0] = cond; }
    if (nlhs > 1) { plhs[1] = cut; }
    if (nlhs > 2) { plhs[2] = vol; }
    if (nlhs > 3) { plhs[3] = s; }
    if (nlhs > 4) { plhs[4] = cc; }
    if (nlhs > 5) { plhs[5] = t; }
    mxAssert(nlhs <= 6, "Too many output arguments");
    std::vector< std::pair<size_t, size_t> > clusters;
    triangle_clusters_greedy(mat, cond, cut, vol, s, cc, t, clusters, verbose);

    if (nlhs > 6) { 
        mxArray* cassign = mxCreateDoubleMatrix(clusters.size(),2,mxREAL);
        plhs[6] = cassign;

        double *ci = mxGetPr(cassign);
        double *cv = mxGetPr(cassign) + clusters.size();
        for (size_t i=0; i<clusters.size(); ++i) {
            ci[i] = (double)(clusters[i].first + 1);
            cv[i] = (double)(clusters[i].second + 1);
        }
    }
}