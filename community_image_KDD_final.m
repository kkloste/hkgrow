%% This script makes one final image for KDD 
% based on our preliminary version

% This has moved into the repo now.
load dblp

%%
comnum = 1307;

%%
com = find(C(:,comnum));
vi = 1;
v = com(vi);

[bsethk,condhk,cuthk,volhk] = hkgrow1(A,v,'t',5);
[bsetpr,condpr,cutpr,volpr] = pprgrow(A,v);
