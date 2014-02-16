function core_figure(ginfo,ncp,cc,showtitle,showlbls)
clf;

colors=[166, 206, 227; 31, 120, 180; 178, 223, 138; 51, 160, 44; ]/255;
n = ginfo.nverts;

[h1,h2]=ncpplot(ginfo.neigh.cond,ginfo.neigh.size);
hold on;
        set(h1,'Visible','off');
        set(h2,'Color',0.55*[1,1,1]);

[h1,h2] = ncpplot(ncp.cond, min(ncp.size,n-ncp.size));
set(h1,'Visible','off');
set(h2,'LineWidth',0.5);
set(h2,'Marker','x');
set(h2,'LineStyle','-.');
set(h2,'Color',colors(1,:));

[h1,h2] = ncpplot(ginfo.whiskers.cond, ginfo.whiskers.size);
set(h1,'Visible','off');
set(h2,'LineWidth',1.5);
set(h2,'Marker','none');
set(h2,'LineStyle','--');
set(h2,'Color',colors(3,:));

loglog(cc.size,cc.cond,'o-');
cn1 = cc.cores(cc.cuts.order(cc.ind(1)));
cnlast = cc.cores(cc.cuts.order(cc.ind(end)));
text(cc.size(1)/1.05,cc.cond(1),num2str(cn1),...
    'HorizontalAlignment','right');
text(cc.size(end)/1.05,cc.cond(end),num2str(cnlast),...
    'HorizontalAlignment','right');

if ~showlbls
    xlabel('');
    ylabel('');
end
set_figure_size([3,2.5]);


if showtitle
    title({sprintf('%s, %i verts, %i edges', ...
        ginfo.name, ginfo.nverts, ginfo.nedges), ...
        sprintf('\\kappa=%.3f, Cbar=%.3f', ...
        ginfo.global_cc, ginfo.mean_cc)});
end
set(gca,'XScale','log');
set(gca,'YScale','log');
plot(ginfo.fiedler.size, ginfo.fiedler.cond,'ro','MarkerSize',6);

if min(min(ginfo.neigh.cond),ginfo.fiedler.cond) < 1e-2
    ybot = 4e-5;
    set(gca,'YTick',[1e-4,1e-3,1e-2,1e-1,1]);
else
    ybot = 1e-2;
    set(gca,'YTick',[1e-2,1e-1,1]);
end

ylim([ybot,1]);
xlim([1,1e5]);
set(gca,'XTick',[1,1e1,1e2,1e3,1e4,1e5]);
maxdeg = max(ginfo.degrees);
line(maxdeg*[1,1],ylim,'LineWidth',1,'Color','k');
text(maxdeg*1.05,1.5*ybot,{'max','deg'},'VerticalAlignment','bottom');
if ginfo.nverts/2 < 1e5
    if (ginfo.nverts/2)/maxdeg < 20
        ypos=1.5;
    else
        ypos=1.7;
    end
    line(ginfo.nverts/2*[1,1],ylim,'LineWidth',1,'Color','k');
    text(1.05*(ginfo.nverts/2),ypos*ybot,'$\frac{verts}{2}$',...
        'VerticalAlignment','bottom','interpreter','latex');
end
grid on;
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
greygrid(gca,[0.75,0.75,0.75]);
hold off;