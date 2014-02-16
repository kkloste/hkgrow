function [H,stats]=hypercluster(A,varargin)
% Run a hyperclustering using a command line wrapper a
% This command is not thread safe, it always writes to a file 
% in the directory with the command
% H = hypercluster(A); run a hypercluster on a graph A given by its
% adjacency matrix
% H = hypercluster(A,'key',value, ...)
% specify optional arguments:
%   'maxvol' = maxvol as a total number of edges or fraction of the graph
%   edges
%   'seed' = random seed to use
%   'alpha' = pagerank alpha value
%   'expand' = starting size for each cluster
%   'centrality' = vector of centrality scores
%   'overlap' = the number of clusters a vertex can be in before skipping
%     looking for clusters starting from that vertex
%   'two_core' = 1 to enable the hyperclustering just to run on the
%     two_core of the graph, the piece without tree-like additions
%   'maxcond' 
%   'minsize'
%
% [H,stats] = hypercluster(A,...)
%   also returns statistics on the computed clusters
%   ** it does not return statistics on the extra clusters added for
%   isolated vertices **
%   each row of stats corresponds to the cluster with the same row number.
%     stats(i,1) = size
%     stats(i,2) = volume
%     stats(i,3) = conductance
%     stats(i,4) = computation steps (~ total computation work)
%     stats(i,5) = computation support (~ total vertices examined)
%  
%

% David F. Gleich
% Copyright, 2010

% History
% -------
% :2010-10-25: Added stats output


mydir = fileparts(mfilename('fullpath'));
% parse options
options = struct(varargin{:});
mycmd = fullfile(mydir,'bin','hypercluster');
mygraphfile = fullfile(mydir,'hypercluster-graph.smat');
myoutfile = fullfile(mydir,'hypercluster-graph.hcluster');
mystatfile = fullfile(mydir,'hypercluster-graph.stats');

writeSMAT(mygraphfile,A);

mycmd = [mycmd ' ' mygraphfile];
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
if isfield(options,'overlap') 
    mycmd = [mycmd ' ' sprintf('-overlap %i',options.overlap)];
end
if isfield(options,'minsize') 
    mycmd = [mycmd ' ' sprintf('-minsize %i',options.minsize)];
end
if isfield(options,'maxcond') 
    mycmd = [mycmd ' ' sprintf('-maxcond %i',options.maxcond)];
end

if isfield(options,'centrality')
    mycentralityfile = fullfile(mydir,'hypercluster-graph.centrality');
    dlmwrite(mycentralityfile,options.centrality(:));
    mycmd = [mycmd ' ' sprintf('-centrality %s',mycentralityfile)];
end
if isfield(options,'two_core')
    if options.two_core == 1
        mycmd = [mycmd ' ' sprintf('-two_core')];
    end
end
if isfield(options,'full')
    if options.full == 1
        mycmd = [mycmd ' ' sprintf('-full')];
    end
end
if nargout>1
    mycmd = [mycmd ' -statsfile ' mystatfile];
end
mycmd
status = system(mycmd);
if status==0
    H = readSMAT(myoutfile);
    stats = load(mystatfile);
end