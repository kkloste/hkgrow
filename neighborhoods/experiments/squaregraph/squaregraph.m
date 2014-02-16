%% Test squaring a graph to fix local clusters for tech networks
A = load_graph('as-22july06');

load ../../data/ncpdata.mat
spec = ncpdata('as-22july06').ncp;
n = size(A,1);
%% compute clustering and neighborhood clusters
d = full(sum(A,2));
[cond cut vol cc t] = triangleclusters_mex(A);
nt = sum(t);
nwedges = d.*(d-1)/2;
nw = sum(nwedges);
kappa = nt/nw;
clf;
[h1,h2] = ncpplot(cond,d+1);
[h1,h2] = ncpplot(spec.cond,spec.size);
%% compute clustering for the ball graph with two edges
[cond2 cut2 vol2 size2] = ball_clusters(A,2);
clf;
[h1,h2] = ncpplot(cond,d+1); hold on;
[h1,h2] = ncpplot(cond2,size2);
set(h1,'Color','r');
set(h2,'Color','r');

%% compute clustering for the ball graph with two edges
[cond2 cut2 vol2 size2] = ball_clusters(A,3);
clf;
[h1,h2] = ncpplot(cond,d+1); hold on;
[h1,h2] = ncpplot(cond2,size2);
set(h1,'Color','r');
set(h2,'Color','r');

%% Summary
% Okay, so this isn't the best idea :-)