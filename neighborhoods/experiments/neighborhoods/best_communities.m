%% Show the best communities from a variety of metrics

%% Load data and list it by type
load '../../data/gdata.mat';
clear sgdata;
for graph=gdata.keys
    if ~exist('sgdata','var'), sgdata = gdata(graph{1}); 
    else, sgdata(end+1) = gdata(graph{1}); end
end
types = unique({sgdata.type});
tgdata = struct();
for ti=1:numel(types)
    type = types{ti};
    tgraphs = sgdata(strcmp({sgdata.type},type));
    [~,p] = sort([tgraphs.nverts]);
    tdata = [];
    for gi=1:numel(p)
        if gi==1
            tdata = tgraphs(p(gi));
        else
            tdata(end+1) = tgraphs(p(gi));
        end
    end
    tgdata.(type) = tdata;
end

%% Load other precomputed data
load '../../data/ncpdata.mat';
load '../../data/metisdata.mat';

%% Spit out a table
print_row = @(best) ...
    fprintf('%20s & %6.4f & %5i & %6.4f & %5i & %6.4f & %5i & %6.4f & %5i & %6.4f & %5i \\\\\n',...
    best.name, best.neigh.cond, best.neigh.size, ...
    best.fied.cond, best.fied.size, best.ncp.cond, best.ncp.size, ...
    best.whisker.cond, best.whisker.size, best.metis.cond, best.metis.size);
    

% Use the order from the paper to make it easier to paste in
for type={'ca','soc','tech','web','model'}
    type = type{1}; % remove the cell
    tgraphs = {tgdata.(type).name};
    for gi=1:numel(tgraphs)
        name = tgraphs{gi};
        ginfo = gdata(name);
        ncp = ncpdata(name).ncp;
        metis = metisdata(name).metis;
        best = struct();
        best.name = name;
        [best.neigh.cond,ind] = min(ginfo.neigh.cond); 
        best.neigh.size = ginfo.neigh.size(ind);
        
        best.fied.cond = ginfo.fiedler.cond;
        best.fied.size = ginfo.fiedler.size;
        
        [best.ncp.cond,ind] = min(ncp.cond); 
        best.ncp.size = ncp.size(ind);
        
        [best.metis.cond,ind] = min(metis.cond);
        best.metis.size = metis.size(ind);
        
        [best.whisker.cond,ind] = min(ginfo.whiskers.cond);
        best.whisker.size = ginfo.whiskers.size(ind);
        
        print_row(best);
    end
    fprintf('\\midrule \n');
end

