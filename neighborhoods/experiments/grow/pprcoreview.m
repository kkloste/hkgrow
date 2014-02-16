%% Good partitions of fb-A-one?
% the fb-A_one year network seemed to have VERY good partitions,
% from the ppr core communites.  We investigate them here.

A = load_external_graph('traces/anony-interactions-oneyearA-cc');

%%
load '../../data/gdata.mat'

%%
cc = corecuts(A);
n = size(A,1);
%%
ccuts = cc;
nv = numel(ccuts.cond);
for i=nv-11:-1:nv-20
    nset = ccuts.cuts.order(1:ccuts.ind(i));
    if numel(nset) > n/2
        nset = setdiff(1:n,nset);
    end
    [curset,coregrow.cond(i),coregrow.cut(i),coregrow.vol(i)] = ...
        pprgrow(A,nset,'nruns',4,'maxexpand',nnz(A)/3);
    coregrow.size(i) = min(numel(curset),n-numel(curset));
    coregrow.cond(i),coregrow.size(i)
end
%%
% set i=8 seems to indicate what happens
i=8;
nset = ccuts.cuts.order(1:ccuts.ind(i));
if numel(nset) > n/2
    nset = setdiff(1:n,nset);
end
[curset,coregrow.cond(i),coregrow.cut(i),coregrow.vol(i)] = ...
    pprgrow(A,nset,'nruns',4,'maxexpand',nnz(A)/3);
curset = setdiff(1:n,curset);
%%
xy=igraph_draw(A(curset,curset));
%fgplot(A(curset,curset),xy);
