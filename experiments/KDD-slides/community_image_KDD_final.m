%% This script makes one final image for KDD 
% based on our preliminary version

% This has moved into the repo now.
load dblp

%%
cd neighborhoods
addpath(pwd)
cd ..

%% remove small communities
comsize = sum(C);
C = C(:,comsize > 40);

%%
comnum = 3204;

com = find(C(:,comnum));
%vi = 2;
%v = com(vi);
v = 51998

[bsethk,condhk,cuthk,volhk] = hkgrow1(A,v,'t',5);
fprintf('hk:\n');
[numel(intersect(com,bsethk))/numel(com) numel(intersect(com,bsethk))/numel(bsethk)]

[bsetpr,condpr,cutpr,volpr] = pprgrow(A,v,'expands',[100,3000]);
fprintf('pr:\n');
[numel(intersect(com,bsetpr))/numel(com) numel(intersect(com,bsethk))/numel(bsetpr)]


% find the subset of the graph
gsub = union(com,bsethk);
gsub = union(gsub,bsetpr);

assert(numel(gsub) < 2000,'too big!');

Asub = A(gsub,gsub);
xy = igraph_draw(Asub,'kk');

fullxy = zeros(size(A,1),2);
fullxy(gsub,:) = xy;

clf;
gplot(Asub,xy,'.-'); hold on;
plot(fullxy(com,1),fullxy(com,2),'r.','MarkerSize',25);
plot(fullxy(bsethk,1),fullxy(bsethk,2),'g.','MarkerSize',18);
plot(fullxy(bsetpr,1),fullxy(bsetpr,2),'k.','MarkerSize',12);
plot(fullxy(v,1),fullxy(v,1)

%% Final example and demo
% using the full com set

comnum = 3204;
com = find(C(:,comnum));
v = 51998;

[bsethk,condhk,cuthk,volhk] = hkgrow1(A,v,'t',5);
fprintf('hk:\n');
[numel(intersect(com,bsethk))/numel(com) numel(intersect(com,bsethk))/numel(bsethk)]

[bsetpr,condpr,cutpr,volpr] = pprgrow(A,v,'expands',[100,3000]);
fprintf('pr:\n');
[numel(intersect(com,bsetpr))/numel(com) numel(intersect(com,bsethk))/numel(bsetpr)]

gsub = union(com,bsethk);
gsub = union(gsub,bsetpr);

xone = zeros(n,1);
xone(gsub) = 1;
bfs1 = A*xone;

gsub = union(gsub,find(bfs1));

Asub = A(gsub,gsub);
%xy = igraph_draw(Asub,'fr');
%xy = kamada_kawai_spring_layout(Asub,'maxiter',1500);
%xy = kamada_kawai_spring_layout(Asub,'maxiter',1500,'progressive',xy);
xy = fruchterman_reingold_force_directed_layout(Asub,'iterations',1500,'initial_temp',100,'force_pairs','all');
xy = fruchterman_reingold_force_directed_layout(Asub,'iterations',1500,'initial_temp',100,'force_pairs','all','progressive',xy);

fullxy = zeros(size(A,1),2);
fullxy(gsub,:) = xy;

clf;
gplot(Asub,xy,'.-'); hold on;

wh=ones(1,3);
colors={0.85*wh,0.7*wh,0.55*wh,0.4*wh};


cl=colors{1}; lw=16*1.5^(-1);
xym=1.05*fullxy(com,:);is=convhull(xym(:,1),xym(:,2));
h=fill(xym(is,1),xym(is,2),cl);set(h,'LineWidth',lw,'FaceColor',cl,'EdgeColor',cl);

cl=colors{2}; lw=16*1.5^(-1);
xym=1.05*fullxy(bsethk,:);is=convhull(xym(:,1),xym(:,2));
h=fill(xym(is,1),xym(is,2),cl);set(h,'LineWidth',lw,'FaceColor',cl,'EdgeColor',cl);

cl=colors{3}; lw=16*1.5^(-1);
xym=1.05*fullxy(bsetpr,:);is=convhull(xym(:,1),xym(:,2));
h=fill(xym(is,1),xym(is,2),cl);set(h,'LineWidth',lw,'FaceColor',cl,'EdgeColor',cl);



plot(fullxy(com,1),fullxy(com,2),'r.','MarkerSize',25);
plot(fullxy(bsethk,1),fullxy(bsethk,2),'g.','MarkerSize',18);
plot(fullxy(bsetpr,1),fullxy(bsetpr,2),'k.','MarkerSize',12);

%% Save the final xy coords here
xyfinal = xy;

%%
save 'com_image_example_1.mat' Asub gsub bsetpr bsethk com comnum v xyfinal

%%
load com_image_example_1
clf;
xy = xyfinal;
fullxy = zeros(v,2);
fullxy(gsub,:) = xy;
[px,py] = gplot(Asub,xy);
h = plot(px,py,'.-','LineWidth',0.3,'MarkerSize',6,'Color',0.55*[1,1,1]);
axis off;
axis square;
set_figure_size([4,4]);
%print(gcf,'com_image_1_base.eps','-depsc2');
lims = axis;
set(h,'Visible','off'); 
hold on;
h2 = plot(fullxy(com,1),fullxy(com,2),'.','MarkerSize',10,'Color',0.8*[1,0,0]);
axis off;
axis(lims);
%print(gcf,'com_image_1_truecom.eps','-depsc2');
set(h2,'Visible','off');

hold on;
h2 = plot(fullxy(bsetpr,1),fullxy(bsetpr,2),'.','MarkerSize',8,'Color',0.8*[0,1,0]);
axis off;
axis(lims);
%print(gcf,'com_image_1_pr.eps','-depsc2');
set(h2,'Visible','off');

h2 = plot(fullxy(bsethk,1),fullxy(bsethk,2),'.','MarkerSize',8,'Color',0.8*[0,0,1]);
axis off;
axis(lims);
%print(gcf,'com_image_1_hk.eps','-depsc2');


%%
% Final demo
xy = xyfinal;
[px,py] = gplot(Asub,xy);
h = plot(px,py,'.-','LineWidth',0.3,'MarkerSize',6,'Color',0.55*[1,1,1]);
axis off;
axis square;
hold on;
h2 = plot(fullxy(com,1),fullxy(com,2),'.','MarkerSize',10,'Color',0.8*[1,0,0]);
h2 = plot(fullxy(bsetpr,1),fullxy(bsetpr,2),'.','MarkerSize',8,'Color',0.8*[0,1,0]);
h2 = plot(fullxy(bsethk,1),fullxy(bsethk,2),'.','MarkerSize',8,'Color',0.8*[0,0,1]);
plot(fullxy(v,1),fullxy(v,2),'mo','MarkerSize',25);



%% Example 1

comsize = sum(C);
C = C(:,comsize > 10);

comnum = 5;

com = find(C(:,comnum));
vi = 2;