% Produces images of graphs with a low-conductance set
cd ~/Dropbox/hkgrow/code/experiments/KDD-slides
clear; clc; clf;

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
theta = 7*pi/4;
R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
xy = xy*R;

% PRINT ORIGINAL GRAPH
seed = 9;

clf;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.55*[1,1,1]);
hold on;
scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
axis off;
axis square;
set_figure_size([4,4]);
% print(gcf,'graph_plain.eps','-depsc2');


 % compute, refine hk rankings
     [bset,cond,cut,vol, hkvec, np] = hkgrow1(A,seed); disp([vol,cond,cut])
%     hk = [];
%     hk = hkvec;%(sginds');
d = full(sum(A,1)); n = length(d);
Dinv = spdiags( 1./d', 0, n, n); 
P = A*Dinv;
alpha = 0.85;
hk = (speye(n) - alpha*P)\sparse( seed, 1, 1-alpha, n, 1);

hk(seed)= full(sum(hk)/length(hk));
    hk = hk./norm(hk,1); mx = max(hk); hk = hk./mx;
    % try log scale?
    hk = log10(hk+ 1e-8.*ones(length(hk),1));
    hk = hk + ones(length(hk),1).*abs((min(hk)));
    hk = hk./max(hk);
    
    hk = exp(hk); hk = hk./max(hk);
    
    hkmin = min(hk);
    acoeff = 1/(1-hkmin);
    bcoeff = 1 - acoeff;
    hk = acoeff.*hk + bcoeff*ones( length(hk), 1);
    
%     % spread out the values more
%     nnzhk = find(hk>0);
%     hdist = 1-min(hk(nnzhk)); hk = hk./hdist;
%     hk(nnzhk) = hk(nnzhk) - ones(length(nnzhk),1).*(max(hk)-1);
    
% % PRINT rank heatmap
%hk_colors = zeros( length(hk),3);
%hk_colors(:,2) = hk;
hk_colors = hk;
    
% NOW FOR THE ZOOMED IN IMAGES
    % PRINT just the bset and its neighbors
   clf;
   
    inds = find(A(:,bset)*ones(length(bset),1) ); % extract the bset and its neighbors
    B = zeros(size(A,1));
    B(inds,inds) = A(inds,inds);
    [px,py] = gplot(B,xy);
    h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.55*[1,1,1]);

    hold on;


%     % PRINT just the bset with rank heatmap
    %scatter(xy(inds,1),xy(inds,2),10,hk_colors(inds,:),'fill');
    scatter(xy(inds,1),xy(inds,2),10,hk_colors(inds),'fill');
    colormap(winter);
    scatter(xy(seed,1),xy(seed,2),14,[0,0,0],'fill');    
    axis off;
    % axis square;
    hold on;
    
    xlim([-81.5, 35]);
    ylim([35, 118.6]);
    set_figure_size([4,4]);
%    print(gcf,'graphsub_diffusion_aptrank.eps','-depsc2');
    print(gcf,'graphsub_diffusion_aptrank','-dpng','-r600');
    
 