function [error, time_x] = test_gsqexpmseed(filename, path, numnonzeros, tol, cluster_size)
% TEST_GSQEXPMSEED   [error, time] = test_gsqexpmseed(filename, path,numnonzeros,tol,cluster_size)
%   Computes an approximate heatkernel pagerank vector for (I-alpha*P)
%       for preference vector with 'numnonzeros' number of nonzeros
%   chosen uniformly at random.
%   Compares result with explicitly computed hkpr vector,
%   via Higham/Al-Mohy's expmv
%       To test accuracy, the top 'cluster_size' entries
%   are compared with those of the explicitly computed
%   hkpr vector.
%
%   Default values:
%
%       numnonzeros = 20
%       tol = 1e-5
%       cluster_size= 100
%

if (nargin == 0)
    filename = 'email';
    path = '~/data';
end

if (nargin <=2)
numnonzeros = 20;
tol = 1e-5;
cluster_size = 100;
end

if (nargin <=3)
tol = 1e-5;
cluster_size = 100;
end

if (nargin <=4)
cluster_size = 100;
end

tic; A = load_graph(filename,path); P = colnormout(A); clear A; n = size(P,1); prep_time = toc


c = randi(n,numnonzeros,1); ec = zeros(n,1); ec(c)=1;
tic; [hkpr,sdummy,m,mv,mvd] = expmv(1,P,ec,[],'single'); time_hkpr = toc

tic; x = gsqexpmseed_mex(P,c,tol); time_x = toc, error = norm(hkpr-x,1)/norm(hkpr,1)
[xval xper] = sort(x,'descend'); [hkprval hkprper] = sort(hkpr, 'descend'); cc =length(intersect(xper(1:cluster_size),hkprper(1:cluster_size)))