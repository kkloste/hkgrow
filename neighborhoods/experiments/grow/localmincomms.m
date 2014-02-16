%% Show the locally minimal neighborhood communities on the ncp

load '../../data/gdata.mat';
load '../../data/ncpdata.mat';

%%

graphs = {'itdk0304-cc'}; %only one figure
for gi=1:numel(graphs)
    %%
    ginfo = gdata(graphs{gi}); 
    gname=ginfo.name;
    if ginfo.nedges > 5e6, fprintf('Skipping %s\n', ginfo.name); continue; end
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
    else
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
   
    
    localmin_figure(A,ginfo,0,0); 
    print(sprintf('localmin-%s.eps',gname),'-depsc2','-painters');
end


%% Determine the average number of local mins that meet our criteria
graphs = gdata.keys;
results = [];
for gi=1:numel(graphs)
    %%
    ginfo = gdata(graphs{gi}); ginfo.name
    %if ginfo.nedges > 5e6, fprintf('Skipping %s\n', ginfo.name); continue; end
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
    else
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
    
    minverts = neighborhoodmin(A,ginfo.neigh.cond,0);
    minverts = minverts(ginfo.degrees(minverts)>=5);
    
    r = struct();
    r.graph = ginfo.name;
    r.src = ginfo.src;
    r.n = ginfo.nverts;
    r.nminverts = numel(minverts);
    
    if isempty(results), results = r; else results(end+1) = r; end
end
%% Play with results to get an interesting number
Tdata = [results.nminverts; results.n];
Tdata(:,Tdata(2,:)>85000)
mean(Tdata(1,:)./Tdata(2,:))
% So we get a 3% average over the nodes with at least 85000 vertices,
% for networkA, we need only consider 17547 seeds, and for 
% soc-LiveJournal, there are only 96,000 seeds. 

%% What about increasing the min number

ginfo = gdata('soc-LiveJournal1'); ginfo.name
%if ginfo.nedges > 5e6, fprintf('Skipping %s\n', ginfo.name); continue; end
if strcmp(ginfo.src,'local')
    A = load_graph(graphs{gi});
else
    A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
end

minverts = neighborhoodmin(A,ginfo.neigh.cond,0);
minverts = minverts(ginfo.degrees(minverts)>=9);

[numel(minverts), size(A,1)]

% So half of the seeds for livejournal are eliminated when
% slightly increasing the minsize to 9


