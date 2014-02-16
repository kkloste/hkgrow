function [H,stats]=spectralncp(A,varargin)
% Run local pagerank clusters using a command line wrapper a C++ code.
%
% This command is not thread safe, it always writes to a file 
% in the directory with the command
%
% H = spectralncp(A); run repeated PageRank clusterings 
% on a graph A given by its adjacency matrix
%
% H = spectralncp(A,'key',value, ...)
% specify optional arguments:
%
%   'maxvol' = maxvol as a total number of edges
%      or fraction of the graph edges
%   'alpha' = pagerank alpha value
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

% David F. Gleich
% Copyright, 2011

% History
% -------
% 2011-10-13: Initial coding

n = size(A,1);

% parse options
options = struct(varargin{:});

pdata.mydir = fileparts(mfilename('fullpath'));
mycmd = fullfile(pdata.mydir,'bin','hypercluster');
mygraphfile = fullfile(pdata.mydir,'spectralncp-graph.smat');
pdata.myoutfile = fullfile(pdata.mydir,'spectralncp-graph.hcluster');
pdata.mystatfile = fullfile(pdata.mydir,'spectralncp-graph.stats');

writeSMAT(mygraphfile,A);

mycmd = [mycmd ' ' mygraphfile ' -method multipagerank'];
mycmd = [mycmd ' -expandfactor 1'];
if isfield(options,'maxvol') 
    mycmd = [mycmd ' ' sprintf('-maxvol %f',options.maxvol)];
else
    mycmd = [mycmd ' ' sprintf('-maxvol 0.50')];
end
if isfield(options,'alpha') 
    mycmd = [mycmd ' ' sprintf('-alpha %f',options.alpha)];
end
if isfield(options,'minsize') 
    mycmd = [mycmd ' ' sprintf('-minsize %i',options.minsize)];
end
if isfield(options,'maxcond') 
    mycmd = [mycmd ' ' sprintf('-maxcond %i',options.maxcond)];
else
    mycmd = [mycmd ' ' sprintf('-maxcond 0.95')];
end
if isfield(options,'maxvisits') 
    mycmd = [mycmd ' ' sprintf('-maxvisits %i',options.maxvisits)];
end

if nargout>1
    mycmd = [mycmd ' -statsfile ' pdata.mystatfile];
end


% mycmd now has the prototype command
% run this for a few values of expand

% run this for extremely local clusters first
expand = 2;
if nargout>1
    [H,stats] = run_multicluster(mycmd,pdata,expand,1:n);
else
    H = run_multicluster(mycmd,pdata,expand,1:n);
end

d = full(sum(A,2));
meand = mean(d);
nedges = nnz(A);

seq = [3 4 5 10 15 20];
epsvals = [seq 10*seq 100*seq 1000*seq];

for epsval = epsvals
    verts = randperm(n);
    if epsval*meand^2 > nedges
        % reduce the number of seeds
        verts = verts(1:min(ceil(0.1*n),1e5));
        if epsval*meand > nedges
            verts = verts(1:min(ceil(0.1*length(verts)),1e4));
        end
    end
    
    lastrun = false;
    if all(epsval*d(verts) >= nedges)
        lastrun = true;
    end
    
    if nargout>1
        [H1,stats1] = run_multicluster(mycmd,pdata,epsval,verts);
        stats = [stats; stats1];
    else
        H1 = run_multicluster(mycmd,pdata,epsval,verts);
    end
    H = [H; H1];

    if lastrun
        break
    end
end

function [H,stats]=run_multicluster(mybasecmd,pdata,expand,vertices)

mycmd = mybasecmd;

mycmd = [mycmd ' ' sprintf('-expand %f',expand)];

myvertexfile = fullfile(pdata.mydir,'spectralncp-graph.vertices');
dlmwrite(myvertexfile,vertices(:)-1,'precision','%i');
mycmd = [mycmd ' ' sprintf('-vertices %s',myvertexfile)];


status = system(mycmd);
if status==0
    H = readSMAT(pdata.myoutfile);
    if nargout>1
      stats = load(pdata.mystatfile);
    end
end
