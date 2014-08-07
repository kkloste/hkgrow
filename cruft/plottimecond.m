% first load runtimes

load timecond;
[vals perm] = sort(inputsize,'ascend');
sortedinputsize = inputsize(perm);
sortedpercdata = percdata(perm,:,:);
hold all
% do hkgrow
fnid = 1;
plot(log10(sortedinputsize),sortedpercdata(:,2,fnid),'g.-');
plot(log10(sortedinputsize),sortedpercdata(:,1,fnid),'g.--');
plot(log10(sortedinputsize),sortedpercdata(:,3,fnid),'g.--');


% do pprgrow
fnid = 2;
plot(log10(sortedinputsize),sortedpercdata(:,2,fnid),'b.-');
plot(log10(sortedinputsize),sortedpercdata(:,1,fnid),'b.--');
plot(log10(sortedinputsize),sortedpercdata(:,3,fnid),'b.--');


title('Conductance*time: hk vs. ppr');
xlabel('log10(|V|+|E|)');
ylabel('Conductance*time');
legend('hkgrow 50%', '25%', '75%', 'pprgrow 50%', '25%', '75%');
set_figure_size([5,3]);
print(gcf,strcat('runtimes','.eps'),'-depsc2');

