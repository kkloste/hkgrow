randgauss = @(sigma,offset) @(n,j) sigma*randn(n,2)+repmat(offset,n,1);
g1 = randgauss(1,[0,0]);
g2 = randgauss(0.5,[-1,-1]);
g3 = randgauss(0.25,[1,-1]);
randfunc = @(n,j) [g1(n/4,j); g2(n/2,j); g3(n/4,j)];
xy = randfunc(1000,2);
plot(xy(:,1),xy(:,2),'.')'
[A,xy] = randgeom(1000, randfunc);
gplot(A,xy);
xyorig = xy;
%%
%xy = igraph_draw(A,'kk');
%[V,D] = eigs(laplacian(A),3,'SA');
%xy = V(:,[2:3]);
xy = xyorig;
gplot(A,xy);

%%
save 'rgeom-ex1.mat' A xy seed;

%%
seed = 76;
t1 = 3;
t2 = 15;
alpha1 = 0.85;
alpha2 = 0.99;
tols = [1e-2,1e-3,1e-4,1e-5];
expands = [100,1000,10000,100000];

clf;
for i=1:2
    switch i
        case 1, alpha = alpha1; t = t1;
        case 2, alpha = alpha2; t = t2;
    end
    for j = 1:numel(tols)
        e = expands(j);
        tol = tols(j);
        subplot3(4,4,2*(i-1)+1,j,0)
        [S,prcond] = pprgrow1(A,seed,'alpha',alpha,'expand',e);
        gplot(A,xy,'.-');
        hold on; scatter(xy(S,1),xy(S,2),10,'r','filled'); hold off
        axis off;
        title(sprintf('cond = %f\n(%.2f,%.0e)', prcond,alpha,1./e));

        subplot3(4,4,2*(i-1)+2,j,0)
        [S,hkcond] = hkgrow1(A,seed,'t',t,'eps',tol);
        gplot(A,xy,'.-');
        hold on; scatter(xy(S,1),xy(S,2),10,'r','filled'); hold off
        title(sprintf('cond = %f\n(%i,%.0e)', hkcond,t,tol));
        axis off;
        drawnow;
    end
end

%% Final picture split into two.
seed = 76;
t1 = 3;
t2 = 15;
alpha1 = 0.85;
alpha2 = 0.99;
tols = [1e-2,1e-3,1e-4,1e-5];
expands = [100,1000,10000,100000];
tx = (max(xy(:,1)) - min(xy(:,1)))/2+min(xy(:,1));
ty = max(xy(:,2))+0.65;
for i=1:2
    switch i
        case 1, alpha = alpha1; t = t1;
        case 2, alpha = alpha2; t = t2;
    end
    figure(i); clf;
    for j = 1:numel(tols)
        e = expands(j);
        tol = tols(j);
        subplot3(2,4,1,j,-0.025,0,0)
        [S,prcond] = pprgrow1(A,seed,'alpha',alpha,'expand',e);
        [px,py] = gplot(A,xy);
        plot(px,py,'k.-','LineWidth',0.1,'MarkerSize',7)
        hold on; scatter(xy(S,1),xy(S,2),10,'r','filled'); hold off
        axis off;
        text(tx,ty,sprintf('cond = %f', prcond),'HorizontalAlignment','center');

        subplot3(2,4,2,j,-0.025,0,0)
        [S,hkcond] = hkgrow1(A,seed,'t',t,'eps',tol);
        [px,py] = gplot(A,xy);
        plot(px,py,'k.-','LineWidth',0.1,'MarkerSize',7)
        hold on; scatter(xy(S,1),xy(S,2),10,'r','filled'); hold off
        text(tx,ty,sprintf('cond = %f', hkcond),'HorizontalAlignment','center');
        axis off;
        drawnow;
    end
    set_figure_size([8,3.25]);
    print(gcf,sprintf('geomfig-%i.eps',i),'-depsc2','-painters');
end

%%
figure(3);
plot(px,py,'k.-','LineWidth',0.25,'MarkerSize',10)
plot(xy(seed,1),xy(seed,2),'ro','MarkerSize',20);
set_figure_size([6,6]);
axis off;
print(gcf,'geomgraph.eps','-depsc2');


