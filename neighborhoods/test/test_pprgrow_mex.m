A = load_graph('dolphins');
[H,stats] = multicluster(A,1:size(A,1),...
    'expandfactor',1,'expand',10,'minsize',1,'maxvol',1,'maxcond',1,'alpha',0.99);
d = full(sum(A,2));
for i=1:size(A,1)
    % multicluster uses (d(i)+1)*expand to set the target volume
    bestclus = pprgrow_mex(A,i,(d(i)+1)*10,0.99);
    assert(isempty(setdiff(bestclus,find(H(i,:)))));
    assert(isempty(setdiff(find(H(i,:)),bestclus)));
end