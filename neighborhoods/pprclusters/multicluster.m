function [H,stats]=multicluster(A,vertices,varargin)
% Run a multiclustering using a command line wrapper a C++ code.
%
% This command is not thread safe, it always writes to a file 
% in the directory with the command
%
% H = multicluster(A,vertices); run repeated PageRank clusterings 
% on a graph A given by its adjacency matrix, where each cluster starts
% from a vertex in vertices.
%
% H = multicluster(A,vertices,'key',value, ...)
% specify optional arguments:
%
%   'maxvol' = maxvol as a total number of edges
%      or fraction of the graph edges
%   'seed' = random seed to use
%   'alpha' = pagerank alpha value
%   'expand' = starting size for each cluster
%   'expandfactor' = the random graph factor for the expansion size
%   'maxcond' = the largest conductance cluster to consider
%   'minsize' = the smallest cluster we'll consider
%
% [H,stats] = multicluster(A,...)
%
%   also returns statistics on the computed clusters
%   ** it does not return statistics on the extra clusters added for
%   isolated vertices **
%   each row of stats corresponds to the cluster with the same row number.
%     stats(i,3) = size
%     stats(i,4) = volume
%     stats(i,5) = conductance
%     stats(i,6) = computation steps (~ total computation work)
%     stats(i,7) = computation support (~ total vertices examined)
%  
%

% David F. Gleich
% Copyright, 2011

% History
% -------



mydir = fileparts(mfilename('fullpath'));
% parse options
options = struct(varargin{:});
mycmd = fullfile(mydir,'bin','hypercluster');
mygraphfile = fullfile(mydir,'multicluster-graph.smat');
myoutfile = fullfile(mydir,'multicluster-graph.hcluster');
mystatfile = fullfile(mydir,'multicluster-graph.stats');

writeSMAT(mygraphfile,A);

mycmd = [mycmd ' ' mygraphfile ' -method multipagerank'];
if isfield(options,'maxvol') 
    mycmd = [mycmd ' ' sprintf('-maxvol %f',options.maxvol)];
end

if isfield(options,'seed') 
    mycmd = [mycmd ' ' sprintf('-seed %i',options.seed)];
end

if isfield(options,'alpha') 
    mycmd = [mycmd ' ' sprintf('-alpha %f',options.alpha)];
end
if isfield(options,'expand') 
    mycmd = [mycmd ' ' sprintf('-expand %f',options.expand)];
end
if isfield(options,'expandfactor') 
    mycmd = [mycmd ' ' sprintf('-expandfactor %f',options.expandfactor)];
end
if isfield(options,'minsize') 
    mycmd = [mycmd ' ' sprintf('-minsize %i',options.minsize)];
end
if isfield(options,'maxcond') 
    mycmd = [mycmd ' ' sprintf('-maxcond %i',options.maxcond)];
end

myvertexfile = fullfile(mydir,'multicluster-graph.vertices');
dlmwrite(myvertexfile,vertices(:)-1);
mycmd = [mycmd ' ' sprintf('-vertices %s',myvertexfile)];

if nargout>1
    mycmd = [mycmd ' -statsfile ' mystatfile];
end
status = system(mycmd);
if status==0
    H = readSMAT(myoutfile);
    if nargout>1
      stats = load(mystatfile);
    end
end
