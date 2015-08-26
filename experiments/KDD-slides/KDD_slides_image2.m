% Produces images of graphs with a low-conductance set
cd ~/Dropbox/hkgrow/code/experiments/KDD-slides
clear; clc; clf;
%%
load com_image_example_1

addpath ../..; % so hkgrow can be used.


% extract smaller subgraph (more easily visualized)
n = size(Asub,1);
bounds = [-75, 165, 190, -190];
sginds = [];
for j = 1:length(xyfinal),
     x = xyfinal(j,1); y = xyfinal(j,2);
    if x < bounds(2) && x > bounds(1) && y < bounds(3) && y > bounds(4),
        sginds(end+1)=j;
    end
end

A=[];
A = Asub(sginds,sginds);
xy = xyfinal(sginds,:);

% rotation matrix, so graph looks better
theta = pi/2;
R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
xy = xy*R*R*R;


%% PRINT ORIGINAL GRAPH
seed = 9;

clf;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.55*[1,1,1]);
hold on;
scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
axis off;
axis square;
set_figure_size([4,4]);
 print(gcf,'graph_plain.eps','-depsc2');


 % compute, refine hk rankings
    [bset,cond,cut,vol, hkvec, np] = hkgrow1(A,seed); disp([vol,cond,cut])
    hk = [];
    hk = hkvec;%(sginds');
    hk = hk./norm(hk,1); mx = max(hk); hk = hk./mx;
    % try log scale?
    hk = log10(hk+ 1e-8.*ones(length(hk),1));
    hk = hk + ones(length(hk),1).*abs((min(hk)));
    hk = hk./max(hk);
    % spread out the values more
    nnzhk = find(hk>0);
    hdist = 1-min(hk(nnzhk)); hk = hk./hdist;
    hk(nnzhk) = hk(nnzhk) - ones(nnzhk,1).*(max(hk)-1);
    
% PRINT rank heatmap
    scatter(xy(:,1),xy(:,2),10,hk,'fill');
    colormap(winter);
    scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
    axis off;
    axis square;
    set_figure_size([4,4]);
    hold on;
    print(gcf,'graph_diffusion.eps','-depsc2');
    
% PRINT bset
    scatter(xy(bset,1),xy(bset,2),10,[1,0.8,0],'fill');
    scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
    set_figure_size([4,4]);
    hold on;
    print(gcf,'graph_bset.eps','-depsc2');
    

% PRINT SEED's NEIGHBORHOOD
clf;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.55*[1,1,1]);
hold on
    scatter(xy(find(A(seed,:)),1),xy(find(A(seed,:)),2),10,[0.8,0.3,0.3],'fill');
    scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
axis off;
axis square;
set_figure_size([4,4]);
% hold on;
 print(gcf,'graph_seed_neighb.eps','-depsc2');
 
    
% NOW FOR THE ZOOMED IN IMAGES
    % PRINT just the bset and its neighbors
   clf;
   
    inds = find(A(:,bset)*ones(length(bset),1) ); % extract the bset and its neighbors
    B = zeros(size(A,1));
    B(inds,inds) = A(inds,inds);
    [px,py] = gplot(B,xy);
    h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.55*[1,1,1]);
    axis off;
    axis square;
    set_figure_size([4,4]);
    hold on;
    print(gcf,'graphsub_plain.eps','-depsc2');

%     % PRINT just the bset with rank heatmap
    scatter(xy(inds,1),xy(inds,2),10,hk(inds),'fill');
    colormap(winter);
    scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
    axis off;
    axis square;
    set_figure_size([4,4]);
    hold on;
    print(gcf,'graphsub_diffusion.eps','-depsc2');
    
 %   % PRINT bset
    scatter(xy(bset,1),xy(bset,2),10,[1,0.8,0],'fill');
    scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
    set_figure_size([4,4]);
    hold on;
    print(gcf,'graphsub_bset.eps','-depsc2');
    
    
%
clf;
    inds = find(A(:,bset)*ones(length(bset),1) ); % extract the bset and its neighbors
    B = zeros(size(A,1));
    B(inds,inds) = A(inds,inds);
    [px,py] = gplot(B,xy);
    h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.55*[1,1,1]);
    axis off;
    axis square;
    set_figure_size([4,4]);
    hold on;
    print(gcf,'graphsub_plain.eps','-depsc2');
    
scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
scatter(xy(find(A(seed,:)),1),xy(find(A(seed,:)),2),10,[0.8,0.3,0.3],'fill');
axis off;
axis square;
set_figure_size([4,4]);
% hold on;
 print(gcf,'graphsub_seed_neighb.eps','-depsc2');
