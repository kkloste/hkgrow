%% Covering with good neighborhood communities

load '../../data/gdata.mat';

%%
graphs = gdata.keys;
for gi=1:numel(graphs)
    %%
    ginfo = gdata(graphs{gi}); ginfo.name
    %if ginfo.nedges > 5e6, fprintf('Skipping %s\n', ginfo.name); continue; end
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
    else
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
    coverage_figure(ginfo,A,1) ;
    pause;
end

%% Good examples
% ca-AstroPh
% cond-mat-2005-fix-cc (like dblp-cc)

%% Bad examples
% Penn94
% as-22july06

%% to show in paper
% ca-AstroPh/as-22july06

%% Print all for website
% Remove the title.
graphs = gdata.keys;
for gi=1:numel(graphs)
    %%
    ginfo = gdata(graphs{gi}); 
    gname = ginfo.name
    %if ginfo.nedges > 5e6, fprintf('Skipping %s\n', ginfo.name); continue; end
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
    else
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
    coverage_figure(ginfo,A,0) ;
    print(sprintf('cover-%s.eps',gname),'-depsc2','-painters');
end

