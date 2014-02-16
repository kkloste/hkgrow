%% Show the locally minimal neighborhood communities on the ncp

load '../../data/gdata.mat';
load '../../data/ncpdata.mat';

%%
% In the paper = email-Enron-cc, ca-AstroPh-cc, arxiv,networkA-anonymized
%graphs = {'itdk0304-cc'}; %only one figure
%graphs = {'cond-mat-2005-fix-cc','ca-AstroPh-cc','email-Enron-cc'}; %only one figure
%graphs = {'arxiv-ubc'};
allgraphs= gdata.keys;
graphs = {'email-Enron-cc','ca-AstroPh-cc','arxiv-ubc','rand-hyper-4000',...
    'rand-ff-25000-0.4','networkA-anonymized',allgraphs{:}};
% full set
% web
for gi=1:numel(graphs)
    %%
    ginfo = gdata(graphs{gi}); 
    gname=ginfo.name;
    %if ginfo.nedges > 5e6, fprintf('Skipping %s\n', ginfo.name); continue; end
    if strcmp(ginfo.src,'local')
        A = load_graph(graphs{gi});
    else
        A = load_external_graph(fullfile(ginfo.src,graphs{gi}));
    end
    
    ncp=ncpdata(graphs{gi}).ncp;
    gname = graphs{gi};
    
    cc = corecuts(A);
    pprgrow_figure(A,ginfo,ncp,cc,0,0);
   print(sprintf('pprgrow-%s.eps',gname),'-depsc2','-painters');
end
%%

