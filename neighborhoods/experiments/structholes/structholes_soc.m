A = load_external_graph('traces/anony-interactions-oneyearA-cc');
%%
% Extract the largest neighborhood community
d = full(sum(A,2));
[ignore p] = sort(d,'descend');
ind = p(6);

Nlarg = [ind; find(A(:,ind))];
A1 = A(Nlarg,Nlarg);
xy=igraph_draw(A1,'lgl');
fgplot(A1,xy)


%% 