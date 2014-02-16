%% Test Sesh's idea for getting reasonable clusters

%%
addpath('~/dev/matlab-bgl');

%% Let's look at a co-authorship network
A = load_graph('cond-mat-2005-fix-cc');
assert(all(nonzeros(A)==1));

%% 
d = sum(A,2);
ccs = clustering_coefficients(A);
plot(d,ccs,'.');

%% 
% get the fiedler vector and best normalized cut
[nf,fo,fcuts] = nfiedler(A);
[bestcut,ind] = min(fcuts.conductance);
bestset = fcuts.order(1:ind);
plot(fcuts.conductance);
title(sprintf('bestcut = %f, size=%i',bestcut,min(size(A,1)-ind,ind)));


%% 
% Show the cuts associated with the largest degree vertices
d = sum(A,2);
[~,dorder] = sort(d,'descend');
volG = sum(d);
cond = [];
for i=1:min(100,numel(dorder))
    vert = dorder(i);
    set = [vert; find(A(:,vert))];
    vol = sum(d(set));
    cond(i) = cutsize(A,set)/min(vol,volG-vol);
end
plot(cond);


%%
% These don't seem to be too great.  Can we make them better if we make
% them bigger?
d = full(sum(A,2));
[~,dorder] = sort(d,'descend');
volG = sum(d);
for i=1:10
    bfsd = bfs(A,dorder(i));
    [ignore bfsp] = sort(bfsd);
    cuts = cutsweep(A,bfsp);
    plot(cuts.conductance);
    title(sprintf('vertex %i, degree %i', dorder(i), d(dorder(i))));
    pause(2);
end

%% 
% Consider points based on maximum clustering/degree.  We normalize
% degrees (nd) to be between 0 and 1, where 1 = maxdegree
% then we sort by distance in the cc/nd space to the point (1,1);
d = full(sum(A,2));
ccs = clustering_coefficients(A);
nd = d/max(d);
dccdist = (nd-1).^2 + (ccs-1).^2;
[~,dorder] = sort(dccdist);
volG = sum(d);
cond = [];
for i=1:min(100,numel(dorder))
    vert = dorder(i);
    set = [vert; find(A(:,vert))];
    vol = sum(d(set));
    cond(i) = cutsize(A,set)/min(vol,volG-vol);
end
plot(cond);

%%
% Try just clustering coefficient order.
ccs = clustering_coefficients(A);
d = sum(A,2);
[~,dorder] = sort(ccs,'descend');
volG = sum(d);
cond = [];
for i=1:min(1,numel(dorder))
    vert = dorder(i);
    set = find(A(:,vert));
    vol = sum(d(set));
    cond(i) = cutsize(A,set)/min(vol,volG-vol);
end
plot(cond);

%% Check cores
[cns dorder] = core_numbers(A);
dorder = dorder(end:-1:1);
d = full(sum(A,2));
volG = sum(d);
cond = [];
for i=1:min(100,numel(dorder))
    vert = dorder(i);
    set = find(A(:,vert));
    vol = sum(d(set));
    cond(i) = cutsize(A,set)/min(vol,volG-vol);
end
plot(cond);

