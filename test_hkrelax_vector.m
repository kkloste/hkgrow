%% This function is 

% mex -O -largeArrayDims hktest_mex.cpp
load four_clusters
%hkgrow_mex(A,set,t,eps,debugflag)
hktest_mex(A,1,5,0.01,1)
%%
mex -O -largeArrayDims hkseed_mex.cpp
[set,cond,cut,vol,x] = hkseed_mex(A,2,5,0.01,1);
