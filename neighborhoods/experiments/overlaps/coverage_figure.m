function coverage_figure(ginfo,A,showtitle)
clf;
colors=[166, 206, 227; 31, 120, 180; 178, 223, 138; 51, 160, 44; ]/255;

cond = ginfo.neigh.cond;
[ignore ppure] = sort(cond);
[ignore pdeg] = sort(cond./ginfo.degrees);
[pure_maxcond,pure_coverage] = neighcover(A,cond,ppure);
[deg_maxcond,deg_coverage] = neighcover(A,cond,pdeg);

loglog(pure_coverage,pure_maxcond,'.-','LineWidth',1.5,'Color',colors(2,:));
hold on;
%loglog(deg_coverage,deg_maxcond,'.-','LineWidth',1.5,'Color',colors(4,:));

xlim([1e-3,1]);
ylim([5e-3,1]);

%legend('Greedy','Deg-Greedy','Location','Southeast');
%legend boxoff;

if showtitle
    title({sprintf('%s, %i verts, %i edges', ...
        ginfo.name, ginfo.nverts, ginfo.nedges), ...
        sprintf('\\kappa=%.3f, Cbar=%.3f', ...
        ginfo.global_cc, ginfo.mean_cc)});
end

% Mark the threshold at 0.1 cond and 0.25 cond
xl = xlim();
threshs = [0.1 0.25 prctile(cond,25)];
covers = {{pure_maxcond,pure_coverage}};
for ti=1:numel(threshs)
    thresh = threshs(ti);
    for ci=1:numel(covers)
        maxcond = covers{ci}{1};
        coverage = covers{ci}{2};
        
        ind = find(maxcond<thresh,1,'last');
        if isempty(ind), continue; end;
        text(coverage(ind),maxcond(ind),sprintf('%.1f%%',coverage(ind)*100),...
            'HorizontalAlignment','left','VerticalAlignment','top');
    end
    text(xl(1)*1.05,thresh,sprintf('\\phi=%.2f',thresh),...
        'HorizontalAlignment','left','VerticalAlignment','bottom');
    line(xl,thresh*[1,1],'LineWidth',1,'Color','k');
end


box off;

grid on;
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
greygrid(gca,[0.75,0.75,0.75]);
hold off;

set_figure_size([3,2.5]);


