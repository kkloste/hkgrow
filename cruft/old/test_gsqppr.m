function [error, time_x] = test_gsqppr(filename, path, numnonzeros, alpha, cluster_size)
% TEST_GSQPPR   [error, time] = test_gsqppr(filename, path)
%   Computes an approximate pagerank vector for (I-alpha*P)
%       for preference vector with 'numnonzeros' number of nonzeros
%   chosen uniformly at random.
%   Compares result with explicitly computed pr vector.
%       To test accuracy, the top 'cluster_size' entries
%   are compared with those of the explicitly computed
%   pr vector.
%
%   Default values:
%
%       numnonzeros = 20
%       alpha       = .9
%       cluster_size= 100
%

if (nargin == 0)
filename = 'email';
path = '~/data';
end


if (nargin <=2)
numnonzeros = 20;
alpha = .9;
cluster_size = 100;
end

if (nargin ==3)
alpha = .9;
cluster_size = 100;
end

if (nargin ==4)
cluster_size = 100;
end

tic; A = load_graph(filename,path); P = colnormout(A); n = size(P,1); A = (sparse(eye(n)).*(1/alpha)-P); prep_time = toc


c = randi(n,numnonzeros,1); ec = zeros(n,1); ec(c)=1; tic; pr = A\(ec.*(1/alpha)); time_pr = toc
tic; x = gsqppr_mex(P,c,alpha,power(10,-5)); time_x = toc, error = norm(pr-x,1)/norm(pr,1)
[xval xper] = sort(x,'descend'); [prval prper] = sort(pr, 'descend'); cc =length(intersect(xper(1:cluster_size),prper(1:cluster_size)))