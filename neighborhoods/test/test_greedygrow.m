A = load_graph('dolphins');
A = A - diag(diag(A));
Gvol = full(sum(sum(A)));
n = size(A,1);
d = sum(A,2) - diag(A);
[cond cut vol s cc t cdata] = triangleclustersgreedy_mex(A);

for i=1:n
    ci = [find(A(:,i)); i];
    verts = greedyclustergrow(A,ci,Gvol,d(i));
    verts_true = cdata(cdata(:,1)==i,2);
    assert(numel(setdiff(verts_true,verts))==0);
    assert(s(i) == length(verts));
    fprintf('dolphins: vertex %i okay!\n', i);
end

%%
A = load_graph('lesmis');
A = A - diag(diag(A));
Gvol = full(sum(sum(A)));
n = size(A,1);
d = sum(A,2) - diag(A);
[cond cut vol s cc t cdata] = triangleclustersgreedy_mex(A);

for i=1:n
    ci = [find(A(:,i)); i];
    verts = greedyclustergrow(A,ci,Gvol,d(i));
    verts_true = cdata(cdata(:,1)==i,2);
    assert(numel(setdiff(verts_true,verts))==0);
    assert(s(i) == length(verts));
    fprintf('lesmis: vertex %i okay!\n', i);
end