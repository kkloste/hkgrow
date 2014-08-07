
% first load heavytwitterp, or pprheavytwitterp



load hktwitterpheavy;

k = 1;

%% size of cluster
scatter([k:length(setsizes)],log10(setsizes(k:length(setsizes))),5,'blue');

hold all;


load pprtwitterpheavy;

scatter([k:length(setsizes)],log10(setsizes(k:length(setsizes))),5,'red');

titlename = strcat('Twitter clustersizes for large-degree seeds');
title(titlename);
xlabel('(1+2000*(x-1))th largest degree node');
ylabel('log10(clustersize)');
legend('hkgrow','pprgrow');

set_figure_size([3,3]);
print(gcf,strcat('twitclus','.eps'),'-depsc2');

