%%
load '../../data/gdata.mat'
load '../../data/ncpdata.mat'

%%
graphs = {'dblp-cc','Penn94'};
%graphs = {'rand-ff-25000-0.4','rand-ff-25000-0.49'};
%graphs = {'itdk0304-cc'};
for gi=1:numel(graphs)
    %%
    ginfo = gdata(graphs{gi});
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
    else
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
    ncp=ncpdata(graphs{gi}).ncp;
    
    %%
    cc = corecuts(A);
    core_figure(ginfo,ncp,cc,1,1);

    %%
    pause;
end

%% Output the graphs
% graphs for paper: graphs = {'dblp-cc','Penn94'};
% but we now output all graphs for the website
graphs = gdata.keys;
for gi=1:numel(graphs)
    ginfo = gdata(graphs{gi});
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
        elsels
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
    ncp=ncpdata(graphs{gi}).ncp;
    gname = graphs{gi};
    
    cc = corecuts(A);
    core_figure(ginfo,ncp,cc,0,0);
    print(sprintf('core-%s.eps',gname),'-depsc2','-painters');
end