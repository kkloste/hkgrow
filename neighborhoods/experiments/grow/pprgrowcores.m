function coregrow = pprcoregrow(A,ccuts)
n = size(A,1);
nv = numel(ccuts.cond);
fprintf('Computing pprgrow clustering ... for %i coresets ... \n', nv);
t0 = tic;
coregrow.cond = zeros(nv,1);
coregrow.cut = zeros(nv,1);
coregrow.vol = zeros(nv,1);
coregrow.size = zeros(nv,1);
for i=1:nv
    nset = ccuts.cuts.order(1:ccuts.ind(i));
    if numel(nset) > n/2
        nset = setdiff(1:n,nset);
    end
    [curset,coregrow.cond(i),coregrow.cut(i),coregrow.vol(i)] = ...
        pprgrow(A,nset,'nruns',4,'maxexpand',nnz(A)/3);
    coregrow.size(i) = min(numel(curset),n-numel(curset));
    stopbar=progressbar(i/nv,3);
    if stopbar, error('stopped'); end
end
coregrow.time = toc(t0);
fprintf('%.1f sec\n', coregrow.time);
best = accumarray(coregrow.size,coregrow.cond,[],@min);
besti = find(best);
loglog(besti,best(besti),'co-');