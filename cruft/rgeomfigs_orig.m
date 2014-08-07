[A,xy] = randgeom(1e4, @(m,n) randn(m,n));
gplot(A,xy,'.-');

%%
[S,prcond] = pprgrow(A,100,'neighborhood',true);
gplot(A,xy,'.-');
hold on; plot(xy(S,1),xy(S,2),'ro'); hold off
title(sprintf('cond = %f', prcond));
%%
[S,hkcond] = hkgrow(A,100,'neighborhood',true);
gplot(A,xy,'.-');
hold on; plot(xy(S,1),xy(S,2),'ro'); hold off
title(sprintf('cond = %f', hkcond));

%%
vert = 100;
neighborhood = true;
save 'rgeom-ex1.mat' A xy vert neighborhood;

%%
[A,xy] = randgeom(1e4, @(m,n) exp(randn(m,n)));
gplot(A,xy,'.-');
%%
[S,prcond] = pprgrow(A,500);
hold on; plot(xy(S,1),xy(S,2),'ro'); hold off
title(sprintf('cond = %f', hkcond));

%%
[A,xy] = randgeom(1e4);
%%
tic, [S,prcond] = pprgrow(A,100); toc
gplot(A,xy,'.-');
hold on; plot(xy(S,1),xy(S,2),'ro'); hold off
title(sprintf('cond = %f', prcond));

%%
tic, [S,hkcond] = hkgrow1(A,100,'t',10,'eps',1e-18); toc
gplot(A,xy,'.-');
hold on; plot(xy(S,1),xy(S,2),'ro'); hold off
title(sprintf('cond = %f', hkcond));
