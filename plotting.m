etype = 1; % 1 - rseed, 2 - hseed, 3 - rhood 4 - hhood
etypelabel = {'randseed', 'heavyseed', 'randhood', 'heavyhood'};
for etype = 1:4

%% conductances
rseedconds = conds(:,:,etype);
boxplot(rseedconds');

titlename = strcat(dataname, ': conductances :', char(etypelabel(etype)));
title(titlename);
xlabel('Graph ID');
ylabel('conductances');
set_figure_size([3,3]);
print(gcf,strcat(dataname,'small','conds',char(etypelabel(etype)),'.eps'),'-depsc2');


%% times
rseedtimes = times(:,:,etype);
boxplot(rseedtimes');

titlename = strcat(dataname, ': runtimes :', char(etypelabel(etype)));
title(titlename);
xlabel('Graph ID');
ylabel('runtimes');
set_figure_size([3,3]);
print(gcf,strcat(dataname,'small','times',char(etypelabel(etype)),'.eps'),'-depsc2');

%% size of cluster
rseedclustsize = log10(setsizes(:,:,etype));
boxplot(rseedclustsize');
titlename = strcat(dataname, ': clustersize :', char(etypelabel(etype)));
title(titlename);
xlabel('Graph ID');
ylabel('log10(clustersize)');
set_figure_size([3,3]);
print(gcf,strcat(dataname,'small','clust',char(etypelabel(etype)),'.eps'),'-depsc2');

end

% %% edge-volume of cluster
% rseedvols = log10(vols(:,:,etype));
% boxplot(rseedvols');
% 
% %% graphdegrees
% dens = gsize(:,2)./gsize(:,1);
% plot(dens')
% 
% %% graphsizes
% semilogy(gsize(:,2));