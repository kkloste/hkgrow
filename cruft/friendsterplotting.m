% first load smalldata or smallppr

load /scratch2/dgleich/kyle/results/friendstertrials;
numtrials = length(conds);
useconds = zeros(numtrials,8);
usetimes = useconds;
useconds(:,[1,3,5,7]) = conds;
usetimes(:,[1,3,5,7]) = times;

load pprfriendstertrials;
useconds(:,[2,4,6,8]) = conds;
usetimes(:,[2,4,6,8]) = times;

%% plot conductances
boxplot(useconds);
title('hk v. ppr on friendster: conductances');
xlabel('Experiment type');
ylabel('Conductances');
set_figure_size([3,2]);
print(gcf,strcat('friendstercomparecond','.eps'),'-depsc2');

%% plot conductances
boxplot(log10(usetimes));
title('hk v. ppr on friendster: times');
xlabel('Experiment type');
ylabel('log10(Times)');
set_figure_size([3,2]);
print(gcf,strcat('friendstercomparetime','.eps'),'-depsc2');
