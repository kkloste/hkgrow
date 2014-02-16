WCtest = sparse([0 1 1 0 0; 1 0 1 0 0; 1 1 0 1 0; 0 0 1 0 1; 0 0 0 1 0]);

[cond,cut,vol,s] = whiskercuts(WCtest);
assert(cond==1/3)
assert(cut==1)
assert(vol==3)
assert(s==2)

%%
A = load_graph('netscience-cc');
d = full(sum(A,2));
Gvol = sum(A(:));
[cond,cut,vol,s,f] = whiskercuts(A);
whiskers = unique(nonzeros(f));
for wi=whiskers(:)'
    w = find(f==wi);
    curcut = cutsize(A,w);
    curvol = sum(d(w));
    assert(cut(wi) == curcut);
    assert(vol(wi) == curvol);
    assert(s(wi) == length(w));
    assert(cond(wi) == curcut/min(Gvol-curvol,curvol));
end
