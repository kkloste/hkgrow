% Produces images of graphs with a low-conductance set
cd ~/Dropbox/hkgrow/code/experiments/KDDmadness-slide
clear; clc; clf;
load com_image_example_1

n = size(Asub,1);
Asub(10,:) = zeros(1,n); Asub(:,10) = zeros(n,1); % node 10 is large degree

seed = 238; %238   239   240   241   242   243   244   245   246   247   248   249   250

 bounds = [-35, 25, 50, -40]; % whole subgraph
% bounds = [-30, 20, 50, 13]; % community
% bounds = [-30, 20, 13, -40]; % complement

inds = find(xyfinal(:,1) < bounds(1) );
for j=1:numel(inds)
    ind = inds(j); Asub(ind,:) = zeros(1,n); Asub(:,ind) = zeros(n,1);
end
inds = find(xyfinal(:,1) > bounds(2) );
for j=1:numel(inds)
    ind = inds(j); Asub(ind,:) = zeros(1,n); Asub(:,ind) = zeros(n,1);
end
inds = find(xyfinal(:,2) > bounds(3) );
for j=1:numel(inds)
    ind = inds(j); Asub(ind,:) = zeros(1,n); Asub(:,ind) = zeros(n,1);
end
inds = find(xyfinal(:,2) < bounds(4) );
for j=1:numel(inds)
    ind = inds(j); Asub(ind,:) = zeros(1,n); Asub(:,ind) = zeros(n,1);
end

sginds = [];
for j = 1:length(xyfinal),
     x = xyfinal(j,1); y = xyfinal(j,2);
    if x < 25 && x > -35 && y < 50 && y > -40,
        sginds(end+1)=j;
    end
end

disp({'total', nnz(Asub(sginds,sginds))})

A = Asub;

% rotation matrix
theta = pi/2;
R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
%R = R';

xy = xyfinal*R;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.45*[1,1,1]);
axis off;
axis square;
set_figure_size([4,4]);
hold on;

print(gcf,'cond_eg.eps','-depsc2');


% 
% % SEED WITH NEIGHBORHOOD
% comm = [];
% comm = find(A(seed,:));
% comm = setdiff(comm, [135]);
% 
% n = size(Asub,1);
% zs = setdiff([1:n], comm);
% A = Asub; A(zs,zs) = zeros(length(zs));
% 
% theta = pi/2; % rotation matrix
% R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
% 
% xy = xyfinal*R;
% [px,py] = gplot(A,xy);
% h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',10,'Color',0.8*[0,1,0]);
% axis off;
% axis square;
% set_figure_size([4,4]);
% hold on;
% 
% print(gcf,'cond_eg_neighb.eps','-depsc2');


% ACTUAL COMMUNITY

comm = [];
for j = 1:length(xyfinal),
     x = xyfinal(j,1); y = xyfinal(j,2);
    if x < 25 && x > -35 && y < 50 && y > 13,
        comm(end+1)=j;
    end
end

comm = setdiff(comm, [135,164,251]);

disp({'inside', nnz(Asub(sginds,comm))})
disp({'outside', nnz(Asub(setdiff(sginds,comm),comm))})

n = size(Asub,1);
zs = setdiff([1:n], comm);
A = Asub;
A(zs,zs) = zeros(length(zs));

% remove outlier, 135, and low-degree nodes
outlier = [135,164,251]; k = length(outlier);
A(outlier,:) = zeros(k,n); A(:,outlier) = zeros(n,k);

% rotation matrix
theta = pi/2;
R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
%R = R';

xy = xyfinal*R;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.8*[1,0,0]);
axis off;
axis square;
set_figure_size([4,4]);
hold on;

print(gcf,'cond_eg_comm.eps','-depsc2');


% JUST THE NEIGHBORHOOD, with community highlighted
clf;
comm = [];
comm = find(Asub(seed,:));
comm = setdiff(comm, [135]);

n = size(Asub,1);
zs = setdiff([1:n], comm);
A = Asub;
A(zs,zs) = zeros(length(zs));

% rotation matrix
theta = pi/2;
R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
%R = R';

xy = xyfinal*R;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.45*[1,1,1]);
axis off;
axis square;
set_figure_size([4,4]);
hold on;

% ACTUAL COMMUNITY

comm = [];
for j = 1:length(xyfinal),
     x = xyfinal(j,1); y = xyfinal(j,2);
    if x < 25 && x > -35 && y < 50 && y > 13,
        comm(end+1)=j;
    end
end

comm = setdiff(comm, [135,164,251]);

n = size(Asub,1);
zs = setdiff([1:n], comm);
A = Asub;
A(zs,zs) = zeros(length(zs));

% remove outlier, 135, and low-degree nodes
outlier = [135,164,251]; k = length(outlier);
A(outlier,:) = zeros(k,n); A(:,outlier) = zeros(n,k);

% rotation matrix
theta = pi/2;
R = eye(2).*cos(theta); R(1,2) = -sin(theta); R(2,1) = sin(theta);
%R = R';

xy = xyfinal*R;
[px,py] = gplot(A,xy);
h = plot(px,py,'.-','LineWidth',0.2,'MarkerSize',6,'Color',0.8*[1,0,0]);
axis off;
axis square;
set_figure_size([4,4]);
hold on;

print(gcf,'cond_eg_closeup.eps','-depsc2');

