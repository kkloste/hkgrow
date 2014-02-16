%% Precompute information about the graphs to make other tasks go faster.

%% load the list of graphs
graphs

%%
metisdata = containers.Map();
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
    Gvol = nnz(A);
    
    fprintf('Computing metis ... \n');
    metis = struct();
    t0 = tic;
    metis = hypercluster_metis(A);
    metis.time = toc(t0);
    
    fprintf('Metis in %.1f sec\n', metis.time);
    
    ginfo.metis = metis;
    
    metisdata(ginfo.name) = ginfo;
    ginfo
    
    save 'metisdata.mat' metisdata
end


