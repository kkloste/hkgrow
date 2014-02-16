function localmin(A,ginfo,showtitle,showlbls)
clf;

minverts = neighborhoodmin(A,ginfo.neigh.cond,0);
minverts = minverts(ginfo.degrees(minverts)>=5);
    

colors=[166, 206, 227; 31, 120, 180; 178, 223, 138; 51, 160, 44; ]/255;
n = ginfo.nverts;

[h1,h2]=ncpplot(ginfo.neigh.cond,ginfo.neigh.size);
hold on;
        set(h1,'Visible','off');
        set(h2,'Color',0.55*[1,1,1]);


loglog(ginfo.neigh.size(minverts),ginfo.neigh.cond(minverts),'r.','MarkerSize',3);
bestmin = accumarray(ginfo.neigh.size(minverts),ginfo.neigh.cond(minverts),[],@min);
bestminsize = find(bestmin);
bestmincond = bestmin(bestminsize);
loglog(bestminsize,bestmincond,'ro');

if ~showlbls
    xlabel('');
    ylabel('');
end
set_figure_size([3,2.5]);


if showtitle
    title({sprintf('%s, %i verts, %i edges', ...
        ginfo.name, ginfo.nverts, ginfo.nedges), ...
        sprintf('\\kappa=%.3f, Cbar=%.3f, %i localmins', ...
        ginfo.global_cc, ginfo.mean_cc, numel(minverts))});
end
set(gca,'XScale','log');
set(gca,'YScale','log');

if min(min(ginfo.neigh.cond)) < 1e-2
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
box off;

grid on;
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
greygrid(gca,[0.75,0.75,0.75]);
hold off;