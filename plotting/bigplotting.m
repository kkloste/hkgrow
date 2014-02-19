% first load friendstertrials or ljournaltrials

%% conductances
rseedconds = conds(:,:);
boxplot(rseedconds);

titlename = strcat(dataname, filename, ': conductances ');
title(titlename);
xlabel('Experiment type');
ylabel('conductances');
% legend('randseed','heavyseed','randhood','heavyhood');
set_figure_size([3,3]);
print(gcf,strcat(dataname, filename, 'conds','.eps'),'-depsc2');


%% times
rseedtimes = times(:,:);
boxplot(rseedtimes);

titlename = strcat(dataname, filename,  ': runtimes ');
title(titlename);
xlabel('Experiment type');
ylabel('runtimes');
% legend('randseed','heavyseed','randhood','heavyhood');
set_figure_size([3,3]);
print(gcf,strcat(dataname, filename, 'times','.eps'),'-depsc2');

%% size of cluster
rseedclustsize = log10(setsizes(:,:));
boxplot(rseedclustsize);
titlename = strcat(dataname,  filename, ': clustersize ');
title(titlename);
xlabel('Experiment type');
ylabel('log10(clustersize)');
% legend('randseed','heavyseed','randhood','heavyhood');
set_figure_size([3,3]);
print(gcf,strcat(dataname, filename, 'clust','.eps'),'-depsc2');


% %% edge-volume of cluster
% rseedvols = log10(vols(:,:));
% boxplot(rseedvols');
          %
          % %% graphdegrees
          % dens = gsize(:,2)./gsize(:,1);
          % plot(dens')
                 %
                 % %% graphsizes
                 % semilogy(gsize(:,2));