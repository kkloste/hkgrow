load rgeom-ex1

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

