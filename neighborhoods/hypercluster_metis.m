function [setdata,H]=hypercluster_metis(A,varargin)
% Run metis to produce a hyperclustering

% David F. Gleich
% Copyright, 2010

% :2011-10-30: Updated for simplecluster project

% parse options
options = struct(varargin{:});

maxvol = 0.5;
if isfield(options,'maxvol') 
    maxvol = options.maxvol;
end

minvol = 10;
if isfield(options,'minvol') 
    minvol = options.minvol;
end

seed = 0;
if isfield(options,'seed') 
    seed = options.seed;
end

nruns = 32;
if isfield(options,'nruns') 
    nruns = options.nruns;
end

if maxvol < 1
    maxvol = nnz(A)*maxvol;
end
if minvol < 1
    minvol = nnz(A)*minvol;
end

maxparts = floor(nnz(A)/minvol);
minparts = ceil(nnz(A)/maxvol);

% TODO handle this case better
% what happens if minvol > nnz(A)?
%maxparts = min(maxparts,minparts*3);
if maxparts<minparts, maxparts=minparts; end

if seed == -1
    rs = RandStream('mt19937ar','seed',sum(100*clock));
else
    rs = RandStream('mt19937ar','seed',seed);
end

nparts = unique(round(logspace(log10(minparts),log10(maxparts),nruns)));
if numel(nparts) < nruns
    % add extra runs at small sizes/many parts
    nextra = nruns - numel(nparts);
    estart = nparts(end-1);
    eend = nparts(end);
    nparts = unique(round([nparts(1:end-2) linspace(estart,eend, nextra+2)]));
end
nruns = numel(nparts);


nsets = sum(nparts);
setdata = struct();
setdata.cond = zeros(nsets,1);
setdata.cut = zeros(nsets,1);
setdata.vol = zeros(nsets,1);
setdata.size = zeros(nsets,1);
si = 1;

H = sparse(size(A,1),0);

for r=1:nruns
    %np = randi(rs,[minparts,maxparts],1);
    np = nparts(r);
    mseed = randi(rs,[0,65536],1);
    p = metismex('PartGraphRecursive',A,np,0,[1 1 1],mseed);
    p = p+1;
    
    
    numclusters = max(p);
    C = sparse(1:size(A,1), p, 1, size(A,1), numclusters);
    valid = true(1,numclusters);
    for c=1:numclusters
        curset = find(C(:,c));
        if numel(curset) == 0 || numel(curset) == size(A,1)
            valid(c) = false; 
            continue; 
        end
        [setdata.cond(si) setdata.cut(si) setdata.vol(si)] = ...
            cutcond_mex(A,curset,nnz(A));
        setdata.size(si) = min(size(A,1) - numel(curset),numel(curset));
        si = si+1;
    end
    if nargout>1
        H = [H C(:,valid)];
    end
    
end
H = H';
if si-1 < numel(setdata.size)
    setdata.size = setdata.size(1:si-1);
    setdata.cond = setdata.cond(1:si-1);
    setdata.vol = setdata.vol(1:si-1);
    setdata.cut = setdata.cut(1:si-1);
end

