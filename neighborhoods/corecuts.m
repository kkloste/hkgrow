function cutdata = corecuts(A)
% CORECUTS Compute cuts based on core removal times
%
% cutdata = corecuts(A)
%
% returns 
%   cutdata.cores: the core numbers for each vertex
%   cutdata.ct: the removal time for each vertex
%   cutdata.cuts: the cut data for the order of the core-time
%   cutdata.bestset: the best set from any cut
%   cutdata.bestset_cond:
%   cutdata.cond the best conductance cut for any core number
%  
%   cutdata.cut 
%   cutdata.vol 
%   cutdata.size 
%     are the cut, volume, and size of the sets with the best conductance
%     score identified in cutdata.cond
%   cutdata.ind(i) is the index of the set of vertices in
%     cutdata.cuts.order for the best cluster at each core
%     i.e. cutdata.cuts.order(1:cutdata.ind(i)) is a set with conductance
%     cutdata.cond(i)
%

n = size(A,1);
cutdata = struct();
[cutdata.cores,cutdata.ct] = core_numbers(A);

assert(min(cutdata.cores) == 1);

[ignore p] = sort(cutdata.ct);
cutdata.cuts = cutsweep(A,p);
[val,ind] = min(cutdata.cuts.conductance);

if ind<((n+1)/2)
    cutdata.bestset = cutdata.cuts.order(1:ind);
else
    cutdata.bestset = cutdata.cuts.order(ind+1:end);
end
cutdata.bestset_cond = val;

ncores = max(cutdata.cores);
cond = Inf*ones(ncores,1);
cut = ones(ncores,1);
vol = ones(ncores,1);
ind = ones(ncores,1);
ssize = ones(ncores,1);

for i=1:n-1
    v = p(i);
    core = cutdata.cores(v);
    if cutdata.cuts.conductance(i) < cond(core)
        cond(core) = cutdata.cuts.conductance(i);
        cut(core) = cutdata.cuts.cut(i);
        vol(core) = cutdata.cuts.vol(i);
        ind(core) = i;
        ssize(core) = min(i,size(A,1)-i);
    end
end

inds = find(isfinite(cond));
cutdata.cond = cond(inds);
cutdata.vol = vol(inds);
cutdata.cut = cut(inds);
cutdata.size = ssize(inds);
cutdata.ind = ind(inds);

assert(min(cutdata.cond) == cutdata.bestset_cond);

