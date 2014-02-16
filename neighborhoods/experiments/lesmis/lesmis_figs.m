%% Show figures of clusters on the lesmis network

A = load_graph('lesmis');
data = ncpneighs(A);
xy = kamada_kawai_spring_layout(A);
xy = fruchterman_reingold_force_directed_layout(A,...
    'progressive',xy,'initial_temp',5,'force_pairs','all');

%% Show the Fiedler cluster
bestset = data.fiedler_set;
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
xyf = xy(data.fiedler_set,:);
hold on; 
Acut = cutgraph(A,data.fiedler_set);
fgplot(Acut,xy,'MarkerSize',0,'Color',[0,0,0],'Alpha',1,'Border',0);
h=scatter(xyf(:,1),xyf(:,2),10,'r','Filled');
hold off;
%title(sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)),'FontSize',8);

%%
set_figure_size([1.75 1.75]);
print('lesmis-fiedler.eps','-depsc2');
fprintf('size=%i, cut=%i, cond=%5.3f\n',numel(bestset),cutsize(A,bestset),cutcond(A,bestset));

%% Show the best neighborhood cluster
[~,bestneigh] = min(data.cond);
bestset = [bestneigh,find(A(bestneigh,:))];
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
hold on;
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',[0,0,0],'Alpha',1,'Border',0);
h=scatter(xys(:,1),xys(:,2),10,'r','Filled'); 
h1=scatter(xy(bestneigh,1),xy(bestneigh,2),15,'Filled');
set(h1,'MarkerFaceColor','r');
set(h1,'MarkerEdgeColor','r');

hold off;
%title(sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));

%%
set_figure_size([1.75 1.75]);
print('lesmis-neigh.eps','-depsc2');
fprintf('size=%i, cut=%i, $\\phi$=%4.2f\n',numel(bestset),cutsize(A,bestset),cutcond(A,bestset));

%% Show the best core cluster
bestset = data.core_set;
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
hold on;
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',[0,0,0],'Alpha',1,'Border',0);
h=scatter(xys(:,1),xys(:,2),10,'r','Filled'); 
hold off;
%title(sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));

%%
set_figure_size([1.75 1.75]);
print('lesmis-core.eps','-depsc2');
fprintf('size=%i, cut=%i, $\\phi$=%4.2f\n',numel(bestset),cutsize(A,bestset),cutcond(A,bestset));

%% Grow the best core cluster
bestset = data.core_set;
bestset = pprgrow(A,bestset);
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
hold on;
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',[0,0,0],'Alpha',1,'Border',0);
h=scatter(xys(:,1),xys(:,2),10,'r','Filled'); 
hold off;
title(sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));

%% Show the best greedy cluster
[~,bestneigh] = min(data.greedy_cond);
Gvol = sum(sum(A)); di = full(sum(A(:,bestneigh)));
bestset = greedyclustergrow(A,[bestneigh,find(A(bestneigh,:))],Gvol,di);
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',25);
hold on;
h=scatter(xys(:,1),xys(:,2),100,'Filled'); 
h1=scatter(xy(bestneigh,1),xy(bestneigh,2),200,'Filled');
set(h1,'MarkerFaceColor','r');
set(h1,'MarkerEdgeColor','r');
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',0.8*[1,0,0],'Alpha',1);
hold off;
title(sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));

%% Grow pagerank clusters
[~,bestneigh] = min(data.cond);
bestset = pprgrow(A,bestneigh);
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',25);
hold on;
h=scatter(xys(:,1),xys(:,2),100,'Filled'); 
h1=scatter(xy(bestneigh,1),xy(bestneigh,2),200,'Filled');
set(h1,'MarkerFaceColor','r');
set(h1,'MarkerEdgeColor','r');
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',0.8*[1,0,0],'Alpha',1);
hold off;
title(sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));

%% Grow pagerank clusters from best greedy vert
[H,stats] = spectralncp(A);
[~,bestclus] = min(stats(:,5));
bestset = find(H(bestclus,:));
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
hold on;
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',[0,0,0],'Alpha',1,'Border',0);
h=scatter(xys(:,1),xys(:,2),10,'r','Filled'); 
bestneigh=stats(bestclus,1);
h1=scatter(xy(bestneigh,1),xy(bestneigh,2),15,'Filled');
set(h1,'MarkerFaceColor','r');
set(h1,'MarkerEdgeColor','r');
hold off;
%sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));
%%
set_figure_size([1.75 1.75]);
print('lesmis-ppr.eps','-depsc2');
fprintf('size=%i, cut=%i, $\\phi$=%4.2f\n',numel(bestset),cutsize(A,bestset),cutcond(A,bestset));

%%  Locally minimal communities
ind = neighborhoodmin(A,data.cond);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
hold on;
h1=scatter(xy(ind,1),xy(ind,2),15,'r','Filled');
hold off;

%%
set_figure_size([1.75 1.75]);
print('lesmis-localmin.eps','-depsc2');

%% Best grown Locally minimal communities
ind = neighborhoodmin(A,data.cond);
addpath('../grow');
grow = pprgrowneigh(A,ind);

[~,bestind] = min(grow.cond);
beststart = ind(bestind);
bestset = pprgrow(A,beststart);
xys = xy(bestset,:);
fgplot(A,xy,'MarkerColor',0.4*[1,1,1],'MarkerSize',5,'Alpha',1,'Border',0);
hold on;
Acut = cutgraph(A,bestset);
fgplot(Acut,xy,'MarkerSize',0,'Color',[0,0,0],'Alpha',1,'Border',0);
h=scatter(xys(:,1),xys(:,2),10,'r','Filled'); 
bestneigh=beststart;
h1=scatter(xy(bestneigh,1),xy(bestneigh,2),15,'Filled');
set(h1,'MarkerFaceColor','r');
set(h1,'MarkerEdgeColor','r');
hold off;
%sprintf('cond=%5.3f size=%i',cutcond(A,bestset),numel(bestset)));
%%
set_figure_size([1.75 1.75]);
print('lesmis-pprgrow.eps','-depsc2');
fprintf('size=%i, cut=%i, $\\phi$=%4.2f\n',numel(bestset),cutsize(A,bestset),cutcond(A,bestset));

