%% Precompute information about the graphs to make other tasks go faster.

%% load the list of graphs
graphs

%%
gdata = containers.Map();
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
    ginfo.nverts = size(A,1);
    ginfo.nedges = nnz(A)/2;
    ginfo.degrees = full(sum(A,2));
    
    fprintf('Computing clustering coeffs ... \n');
    neigh = struct();
    t0 = tic;
    [neigh.cond neigh.cut neigh.vol ginfo.clustercoeffs] = triangleclusters(A);
    neigh.time = toc(t0);
    neigh.size = min(size(A,1) - (ginfo.degrees + 1),(ginfo.degrees + 1));
    
    ginfo.neigh = neigh;
    
    ginfo.ntris = ginfo.clustercoeffs.*(ginfo.degrees).*(ginfo.degrees-1)/2;
    assert(all(abs(round(ginfo.ntris)-ginfo.ntris)<1e-8));
    ginfo.tris = round(ginfo.ntris);
    ginfo.wedges = round((ginfo.degrees).*(ginfo.degrees-1)/2);
    assert(mod(sum(ginfo.tris),3)==0);
    ginfo.ntris = sum(ginfo.tris);
    ginfo.nwedges = sum(ginfo.wedges);
    ginfo.global_cc = ginfo.ntris./ginfo.nwedges;
    ginfo.mean_cc = mean(ginfo.clustercoeffs);
    
    fprintf('Computing Fiedler vector ... \n');
    fied = struct();
    t0 = tic;
    [f,~,cuts,bestset,lam2] = nfiedler(A,1e-7);
    fied.time = toc(t0);
    fied.set = bestset;
    fied.size = min(size(A,1)-numel(bestset),numel(bestset));
    fied.cuts = cuts;
    [fied.cond,fiedind] = min(cuts.conductance);
    fied.cut = cuts.cut(fiedind);
    fied.vol = cuts.vol(fiedind);
    fied.x = f;
    
    ginfo.fiedler = fied;
    ginfo.lam2 = lam2;
    
    ginfo.type = graphlist{i,3};
    ginfo.cores = core_numbers(A);
    
    fprintf('computing whiskers ... \n')
    whiskers = struct();
    
    t0 = tic;
    [whiskers.cond,whiskers.cut,whiskers.vol,whiskers.size] = ...
        bestwhiskercuts(A);
    whiskers.time = toc(t0);
    ginfo.whiskers = whiskers;
    
    gdata(ginfo.name) = ginfo;
    ginfo
    
    save 'gdata.mat' gdata
end


