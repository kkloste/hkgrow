%% Precompute information about the graphs to make other tasks go faster.

%% load the list of graphs
graphs

%%
ncpdata = containers.Map();
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
    
    fprintf('Computing spectralncp ... \n');
    ncp = struct();
    t0 = tic;
    [H stats] = spectralncp(A,'minsize',5,'maxvisits',10);
    ncp.time = toc(t0);
    ncp.cond = stats(:,5);
    ncp.size = stats(:,3);
    ncp.vol = stats(:,4);
    ncp.cut = round(ncp.cond.*min(Gvol-ncp.vol,ncp.vol));
    ncp.start = stats(:,1);
    
    fprintf('Spectral ncp in %.1f sec\n', ncp.time);
    
    ginfo.ncp = ncp;
    
    ncpdata(ginfo.name) = ginfo;
    ginfo
    
    save 'ncpdata.mat' ncpdata
end


