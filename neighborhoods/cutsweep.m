function cutvals = cutsweep(A,p)
% CUTSWEEP All cut values induced by a permutation
%
% cutvals = cutsweep(A,p)
%
% cutvals.cut = cut(S)
% cutvals.vol = vol(S)
% cutvals.edges = edges(S) ( = vol(S) - cut(S))
% cutvals.ncut = cut(S)/vol(S) + cut(S)/vol(notS)
% cutvals.conductance = cut(A,B)/min(vol(A),vol(B))
% cutvals.modularity = 
% cutvals.expansion = 
% cutvals.qcut = cut(S)/min(|S|,n-|S|)
%
% These are all indexed so that cutvals.cut(k) is the cut value where S is
% the set of vertices p(1:k).  Thus, cutvals.cut(end) = 0, because there
% is no other side.

% cutvals.avgcut = cut(S)/|S| + cut(S)/(n-|S|)

% TODO
% Implement weighted cuts

if ~isequal(A,A')
    error('cutsweep:invalidArgument','the input graph A must be symmetric');
end

n = size(A,1);
d = full(sum(spones(A),2));
d = d-diag(A);

cutvals = struct;
cutvals.cut = zeros(n,1);
cutvals.vol = zeros(n,1);
cutvals.order = p;

curset = zeros(n,1); % current set indicator
curcut = 0;
curedges = 0;
curvol = 0;

for i=1:n
    vert = p(i); % add p(i) to the current set
    neighs = find(A(:,vert));
    for j=neighs(:)'
        if j==vert, continue; end % ignore self-loops
        if curset(j)
            % this vertex is in the current set
            curedges = curedges + 2;
            curcut = curcut - 1;
        else
            % this vertex is not in the current set
            curcut = curcut + 1;
        end
    end
    curset(vert) = 1;
    curvol = curvol + d(vert);
    assert(curedges == curvol - curcut);
    
    cutvals.cut(i) = curcut;
    cutvals.vol(i) = curvol;
end

volG = full(cutvals.vol(end));
assert(volG == sum(d));

cutvals.edges = cutvals.vol - cutvals.cut;

cutvals.conductance = cutvals.cut./min(cutvals.vol,volG-cutvals.vol);
cutvals.ncut = cutvals.cut./cutvals.vol + cutvals.cut./(volG-cutvals.vol);
cutvals.modularity = (1/volG) * (cutvals.edges - (1/volG)*cutvals.vol.^2);
cutvals.qcut = cutvals.cut./min(1:n,n-1:-1:0)';

