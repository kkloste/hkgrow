
% first load heavytwitterp, or pprheavytwitterp

filename = 'heavytwitterp';
k = 1;

%% conductances
scatter([k:length(conds)],conds(k:length(conds)),5);

titlename = strcat(dataname, filename, ' : conductances ');
title(titlename);
xlabel('(1+2000*(x-1))th largest degree node');
ylabel('conductances');

set_figure_size([3,3]);
print(gcf,strcat(dataname, filename, 'conds','.eps'),'-depsc2');


%% times
scatter([k:length(times)],log10(times(k:length(times))),5);

titlename = strcat(dataname, filename,  ' : runtimes ');
title(titlename);
xlabel('(1+2000*(x-1))th largest degree node');
ylabel('runtimes');

set_figure_size([3,3]);
print(gcf,strcat(dataname, filename, 'times','.eps'),'-depsc2');

%% size of cluster
scatter([k:length(setsizes)],log10(setsizes(k:length(setsizes))),5);

titlename = strcat(dataname,  filename, ' : clustersize ');
title(titlename);
xlabel('(1+2000*(x-1))th largest degree node');
ylabel('log10(clustersize)');

set_figure_size([3,3]);
print(gcf,strcat(dataname, filename, 'clust','.eps'),'-depsc2');

