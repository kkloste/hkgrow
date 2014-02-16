addpath('~/dev/matlab-bgl');
A = load_graph('dolphins');
n = size(A,1);
d = sum(A,2) - diag(A);
t = ntris(A);
cuts = zeros(n,1); vols=cuts;
for i=1:n
    s = unique([i;find(A(:,i))]);
    cuts(i) = cutsize(A,s);
    vols(i) = sum(d(s));
end
tris = ntris(A);

[cond cut vol cc t] = triangleclusters(A);

assert(all(cut == cuts));
assert(all(vol == vols));
assert(all(t == tris));

%%
[cond cut vol cc t] = triangleclusters_mex(A);

assert(all(cut == cuts));
assert(all(vol == vols));
assert(all(t == tris));
assert(all(cond == cuts./min(sum(sum(A))-vols,vols)));