function grow = pprgrowneigh(A,verts)

n = size(A,1);
nv = numel(verts);

grow.cond = zeros(nv,1);
grow.cut = zeros(nv,1);
grow.vol = zeros(nv,1);
grow.size = zeros(nv,1);
t0 = tic;
for i=1:nv
    nset = [verts(i); find(A(:,verts(i)))];
    [curset,grow.cond(i),grow.cut(i),grow.vol(i)] = ...
        pprgrow(A,nset,'nruns',12,'maxexpand',nnz(A)/5);
    grow.size(i) = min(numel(curset),n-numel(curset));
    stopbar=progressbar(i/nv,3);
    if stopbar, error('stopped'); end
end
grow.time = toc(t0);
fprintf('%.1f sec\n', grow.time);