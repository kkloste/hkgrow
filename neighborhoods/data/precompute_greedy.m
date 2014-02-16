%% Precompute information about the graphs to make other tasks go faster.

%% load the list of graphs
graphs

%%
greedydata = containers.Map();
for i=1:size(graphlist,1)
    ginfo = struct();
    if strcmp(graphlist{i,1},'external')
        [gsrc ginfo.name] = fileparts(graphlist{i,2});
        A = load_external_graph(graphlist{i,2});
        ginfo.src = gsrc;
    else
        A = load_graph(graphlist{i,2});
        ginfo.name = graphlist{i,2};
        ginfo.src = 'local';
    end
    ginfo.name
    assert(all(diag(A)==0));
    assert(max(components(A))==1);
    
    fprintf('Computing neighborhood + greedy clusters ... ');
    greedy = struct();
    t0 = tic;
    [greedy.cond,greedy.cut,greedy.vol,greedy.size] = ...
        triangleclustersgreedy_mex(A, 1);
    greedy.time = toc(t0);
    greedy.size = min(size(A,1) - greedy.size, greedy.size);
    fprintf('%.1f sec\n', greedy.time);
    
    ginfo.greedy = greedy;
    
    greedydata(ginfo.name) = ginfo;
    ginfo
    
    save 'greedydata.mat' greedydata
end


