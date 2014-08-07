function [bestset,bestcond,bestcut,bestvol] = hkgrow(A,vert,varargin)
% HKGROW Grow a cluster around a vertex using a heatkernel-pagerank algorithm
%
% [bestset,cond,cut,vol] = hkgrow(A,vert) computes exp(t*P)*v where
% v is either a single column of the identity or a group of columns of the
% identity and then extract a cluster. The algorithm uses various values of
% t and returns the best conductance cluster among any of them. 
%
% ... hkgrow(A,verts,'key',value,'key',value) specifies optional argument
%
% 'tmax' : the largest value of t to use
% 'nruns' : the number of runs, the default is 32, which will look for clusters 
%    t = 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 60, 80
%    
% 'nruns' : specific the number of heatkernel-pagerank runs.  The default is 32, which will
%   look for clusters with 2, 3, 4, 5, 10, 15, 20, 30, 40, 50, 100, 150, 
%   200, 300, 400, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000,
%   15000,  20000, 30000, 40000, 50000, 100000, 150000, 200000,
%   300000 additional edges for each vertex in the original set of vertices.  
%   This sequence will be continuned until nruns is exhausted.
% 'expands' : a custom sequence of target volumes to consider, e.g. to look
%   for clusters with 12 and 18 edges for each vertex in the original set,
%   then call hkgrow(A,verts,'expands',[12 18]).  
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

% David F. Gleich
% Purdue University, 2011
%
%       heat kernel adaptations with
% Kyle Kloster
% Purdue University, 2014

if nargin<5, debugflag=0; end

p = inputParser;
p.addOptional('nruns',32,@isnumeric);
p.addOptional('alpha',0.99,@isnumeric);
p.addOptional('expands',[]);
p.addOptional('maxexpand',Inf,@isnumeric);
p.parse(varargin{:});

expandseq = [2 3 4 5 10 15];
expands = [];
curmod = 1;
while numel(expands) < p.Results.nruns
    expands = [expands curmod*expandseq];
    curmod = curmod*10;
end
expands = expands(1:p.Results.nruns);

if ~isempty(p.Results.expands)
    expands = p.Results.expands;
end

di=full(max(sum(A(:,vert))));

%expands = unique(round(logspace(log10(3),log10(p.Results.maxexpand),p.Results.nruns)));
bestcond = Inf;
bestset = [];
for ei=1:numel(expands)
    curexpand = expands(ei)*numel(vert)+di;
    if curexpand > p.Results.maxexpand, continue; end

if debugflag==1,    fprintf('hkgrow.m: Called hkgrow_mex on set of size=%i with expand=%i\n', numel(vert), curexpand); end

%    [curset cond cut vol] = hkgrow_mex(A,vert,curexpand,p.Results.alpha, tol);
    [curset cond cut vol] = hkgrow_mex(A,vert,curexpand,alphat, tol, debugflag);

if debugflag==1,    fprintf('hkgrow.m: hkgrow_mex done on set of size=%i with expand=%i\n', numel(vert), curexpand); end

    if cond < bestcond
        bestcond = cond;
        bestset = curset;
        bestvol = vol;
        bestcut = cut;
    end
end