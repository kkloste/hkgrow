%% Load a graph
%A = load_graph('email-Enron-cc');
A = load_graph('as-22july06');
data = ncpneighs(A);

%% Look at the clusters with the neighbormin property
plotncpdata(data);
minverts = neighborhoodmin(A,data.cc./data.cond,1);
hold on;
plot(data.d(minverts)+1,data.cond(minverts),'r*');
hold off;

%% Now greedily grow these clusters
% we've already done this, so just pickup the the data
hold on;
[hpts,hline] = ncpplot(data.greedy_cond(minverts),data.greedy_size(minverts));
set(hpts,'Color','m');
set(hline,'Color','m');
hold off;
%% 
% This isn't a verty good approximation of the green/greedy curve

%% Grow them with PPR instead
Gvol = full(sum(sum(A)));
d = full(sum(A,2));
nmv = numel(minverts);
data.min_cond = zeros(nmv,1);
data.min_size = zeros(nmv,1);
for i=1:nmv
    vert = minverts(i);
    [bestset, cond, cut, vol] = pprgrow(A,vert,'nruns',10);
    data.min_cond(i) = cond;
    data.min_size(i) = length(bestset);
end
%%
hold on;
%plot(data.min_size,data.min_cond,'m.','MarkerSize',12);
[hpts,hline] = ncpplot(data.min_cond,data.min_size);
set(hpts,'Color','y');
set(hline,'Color','y');
hold off;
