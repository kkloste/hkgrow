function neighfig(ginfo,showtitle,showpts,showlbls)
hold on;
[h1,h2] = ncpplot(ginfo.neigh.cond, ginfo.neigh.size);
if ~showpts
    set(h1,'Visible','off');
    set(h2,'LineWidth',0.5);
    set(h2,'Marker','.');
end
if ~showlbls
    xlabel('');
    ylabel('');
end
set(h1,'MarkerSize',3);
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