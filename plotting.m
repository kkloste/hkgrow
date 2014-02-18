etype = 1; % 1 - rseed, 2 - hseed, 3 - rhood 4 - hhood
etypelabel = {'randseed', 'heavyseed', 'randhood', 'heavyhood'};
for etype = 1:4

%% conductances
rseedconds = conds(:,:,etype);
boxplot(rseedconds');

titlename = strcat('conductances of small data:', char(etypelabel(etype)));
title(titlename);
xlabel('Graph ID');
ylabel('conductances of smalldata');
set_figure_size([3,3]);
print(gcf,strcat('conds','small',char(etypelabel(etype)),'.eps'),'-depsc2');


%% times
rseedtimes = times(:,:,etype);
boxplot(rseedtimes');

titlename = strcat('runtimes of small data:', char(etypelabel(etype)));
title(titlename);
xlabel('Graph ID');
ylabel('runtimes of smalldata');
set_figure_size([3,3]);
print(gcf,strcat('times','small',char(etypelabel(etype)),'.eps'),'-depsc2');

%% size of cluster
rseedclustsize = log10(setsizes(:,:,etype));
boxplot(rseedclustsize');
titlename = strcat('clustersize of small data:', char(etypelabel(etype)));
title(titlename);
xlabel('Graph ID');
ylabel('log10 of clustersize');
set_figure_size([3,3]);
print(gcf,strcat('clust','small',char(etypelabel(etype)),'.eps'),'-depsc2');

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