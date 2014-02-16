%% Plots of the properties of the neighborhood communities.
load ../../data/gdata

%% Show them all
% Here, we'll pick out a few interesting figures to put into
% the paper based on properties of all of them.

for gname=gdata.keys
    clf;
    ginfo = gdata(gname{1});
    neighfig(ginfo,1,1)
    pause
end 

%% Best figures for our point
% web-Google
% soc-LiveJournal1
% anony-interactions-oneyearA-cc
% cond-math-2005-fix-cc
% dblp-cc
% hollywood-2009-cc


%% Bad figures for our point
% Penn94
% networkA-anonymized
%
% as-22july06
%
% ca-AstroPh-cc
% email-Enron-cc
% rand-ff-25000-0.49

%% Worst figures for our point
% arxiv-ubc
% itdk0304
% rand-ff-25000-0.4
% p2p-Gnutella25

%%
%print_graphs = {'web-Google','soc-LiveJournal1','anony-interactions-oneyearA-cc',...
%    'itdk0304-cc','arxiv-ubc','ca-AstroPh-cc'};
% This is the full set for the websit now.
print_graphs = gdata.keys
for gi=1:numel(print_graphs)
    gname=print_graphs{gi};
    clf;
    neighfig(gdata(gname),0,0,0);
    print(sprintf('neigh-%s.eps',gname),'-depsc2','-painters');
end
