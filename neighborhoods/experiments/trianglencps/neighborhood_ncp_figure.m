function neighborhood_ncp_figure(graphname)
A = load_external_graph(graphname);
data = ncpneighs(A); 
title({graphname,sprintf('nverts=%i   nedges=%i   maxt=%i   meancc=%f',...
    size(A,1),nnz(A)/2,max(data.t),mean(data.cc))});
[val,ind] = max(data.t);
hold on; loglog(data.d(ind)+1,data.cond(ind),'go'); hold off;
figname = strrep(graphname,'/','-');
print(gcf,[figname '.png'],'-dpng','-r200');

