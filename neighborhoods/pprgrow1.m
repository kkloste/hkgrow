function [bestset,bestcond,bestcut,bestvol] = pprgrow1(A,vert,varargin)
% PPRGROW Grow a cluster around a vertex using a PPR algorithm
%
% [set,cond,cut,vol] = pprgrow1(A,vert) solves a PPR problem and
% extracts a cluster for various cluster volumes ranging from the degree
% of the vertex + 2 to the degree of the vertex + 300000.  In other words, 
% the goal is to find a cluster of up to 300000 additional edges nearby
% the seed vertex.  The returned cluster is the one that has the best
% conductance score among all the clusters found.
%
% [bestset,cond,cut,vol] = pprgrow(A,verts) solves a PPR problem that
% tries to find a set that is 2 to 300000 times the size of the input group
% of vertices.  
%
% ... pprgrow(A,verts,'key',value,'key',value) specifies optional arguments
%
% 'expands' : a custom sequence of target volumes to consider, e.g. to look
%   for clusters with 12 and 18 edges for each vertex in the original set,
%   then call pprgrow(A,verts,'expands',[12 18]).  
% 'maxexpand' : sets a limit on the largest value of the expansion
%   possible, in terms of number of edges.  The default value of this is
%   infinite, but a more reasonable value is nnz(A)/2.  This avoids
%   unnecessary work for larger graphs and can be used to uniformly limit
%   the size of the clusters.
%
% For more detail on these three parameters, please look at how they are
% used in the code.
%
% 'alpha' : the value of alpha to use in pagerank.  The default is 0.99.
%   Using a value of 0.999 is reasonable if your goal is to look for larger
%   clusters.

% warning('pprgrow:badDefaults',...
%    ['this file was edited for fair comparisons in the hkgrow package ' ...
%    'please see the neighborhoods package for better defaults']);

% David F. Gleich
% Purdue University, 2011

p = inputParser;
p.addOptional('alpha',0.99,@isnumeric);
p.addOptional('expand',10000);
p.addOptional('neighborhood',false,@islogical);
p.parse(varargin{:});

if p.Results.neighborhood
    neighs = find(A(:,vert));
    vert = union(vert,neighs);
end

bestcond = Inf;
bestset = [];
curexpand = p.Results.expand;
[curset cond cut vol] = pprgrow_mex(A,vert,curexpand,p.Results.alpha);
if cond < bestcond
    bestcond = cond;
    bestset = curset;
    bestvol = vol;
    bestcut = cut;
end
