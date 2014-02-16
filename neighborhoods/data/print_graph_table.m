%% Print graph table for the paper
load gdata
clear sgdata;
for graph=gdata.keys
    if ~exist('sgdata','var'), sgdata = gdata(graph{1}); 
    else, sgdata(end+1) = gdata(graph{1}); end
end

%% print 
% vertices, edges, mean deg, max deg, $\kappa$, $\bar{C}$
%
print_ginfo = @(g) ...
    fprintf('%20s & %7i & %8i & %5.1f & %5i & %5.3f & %5.3f \\\\ \n', ...
    g.name, g.nverts, g.nedges, mean(g.degrees), max(g.degrees), ...
    g.global_cc, g.mean_cc);

fprintf('\\midrule\n')
types = unique({sgdata.type});
for ti=1:numel(types)
    type = types{ti};
    tgraphs = sgdata(strcmp({sgdata.type},type));
    [~,p] = sort([tgraphs.nverts]);
    for gi=1:numel(p)
        print_ginfo(tgraphs(p(gi)));
    end
    fprintf('\\midrule \n');
end