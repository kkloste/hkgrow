%% Create the pathweights figure for the paper
% The point of this figure is to illustrate the difference between the path
% weighting in PageRank and in the heat-kernel as we vary alpha and t.
% Include time-dependent pagerank in this one!
clear; clc;
clf; hold all;

% Draw two alphas
legtext = {};
alphas = [0.85,0.99];
ts = [1, 5, 15];


for ai=1:numel(alphas)
    alpha = alphas(ai);

    prweight0 = 1-alpha;

    ks = 100;
    prweight = zeros(ks,1);
    prweight(1) = prweight0;
    for k=2:ks
        prweight(k) = alpha*prweight(k-1);
    end
    plot(prweight,'b--','LineWidth',1.2,'Color',[0,0,0.8]);
    legtext{end+1} = sprintf('pr - alpha=%.2f',alpha);
end

for ti=1:numel(ts)
    t = ts(ti);
    hkweight0 = exp(-t);
    hkweight = zeros(ks,1);
    hkweight(1) = hkweight0;
    for k=2:ks
        hkweight(k) = t/(k-1)*hkweight(k-1);
    end
    plot(hkweight,'r-','LineWidth',1.2,'Color',[0.8,0,0]);
    legtext{end+1} = sprintf('hk - t=%g',t);
end



% Time dependent pagerank

tdpr_alpha = 0.85;
tdpr_gamma = 5.0;

egamma = exp(-tdpr_gamma);
tdprweight0 = 1-tdpr_alpha + egamma;

ks = 100;
tdprweight = zeros(ks,1);
tdprweight(1) = tdprweight0;
alphak = 1;
sum_gam = 1;
last_sum_gam_term = 1;
for k=2:ks
    alphak = alphak*tdpr_alpha;
    last_sum_gam_term = last_sum_gam_term*tdpr_gamma/(k-1);
    sum_gam = sum_gam + last_sum_gam_term;
    left_piece = (1-tdpr_alpha)*alphak*(1-egamma*sum_gam);
    right_piece = last_sum_gam_term*alphak*egamma;
    tdprweight(k) = left_piece+right_piece;
end
plot(tdprweight,'.-.','LineWidth',1.2,'Color',[0,0.8,0]);
legtext{end+1} = sprintf('tdpr - alpha=%.2f,  gamma=%.2f', tdpr_alpha, tdpr_gamma);
    

text(10,1e-7,'t=1','VerticalAlignment','bottom','HorizontalAlignment','right');
text(20,1e-7,'t=5','VerticalAlignment','bottom','HorizontalAlignment','right');
text(38,1e-7,'t=15','VerticalAlignment','bottom','HorizontalAlignment','right');
text(87,1e-7,'\alpha=0.85','VerticalAlignment','bottom');
text(87,5e-3,'\alpha=0.99','VerticalAlignment','bottom');
ylabel('Weight');
xlabel('Length');
set(gca,'YScale','log');
set(gca,'LineWidth',0.6);
ylim([1e-8,1]);
%legend(legtext{:});

set_figure_size([4,2]);
print('pathweights_all.png', '-dpng','-r600');