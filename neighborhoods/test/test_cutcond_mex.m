A = load_graph('dolphins');
n = size(A,1);
ntrials = 1000;
Gvol = nnz(A);
assert(all(diag(A)==0));
d = full(sum(A,2));
for ti=1:ntrials
    set = unique(randi(n, randi(n, 1)));
    [cond1,cut1,vol1] = cutcond_mex(A, set, A);
    vol2 = sum(d(set));
    cut2 = cutsize(A,set);
    assert(vol1 == vol2);
    assert(cut1 == cut2);
    cond2 = cut2/min(Gvol-vol2,vol2);
    assert(abs(cut1 - cut2)<1e-8);
end