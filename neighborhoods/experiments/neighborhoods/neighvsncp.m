%% Plots of the properties of the neighborhood communities.
load '../../data/gdata.mat'
load '../../data/ncpdata.mat'

%% Show them all
% Here, we'll pick out a few interesting figures to put into
% the paper based on properties of all of them.

for gname=gdata.keys
    clf;
    ginfo = gdata(gname{1});
    neighncpfig(ginfo,ncpdata(gname{1}).ncp,1,1,1)
    pause
end 

%% Best figures for our point
% email-Enron-cc (bad)
% cond-mat-2005-fix-cc (good)


%%
%print_graphs = {'email-Enron-cc','cond-mat-2005-fix-cc',...
%    'as-22july06','hollywood-2009-cc'};
% That was the initial set og figures for the paper.  THe next
% set is for the website, which is just all of them!
print_graphs = gdata.keys;
for gi=1:numel(print_graphs)
    gname=print_graphs{gi};
    clf;
    neighncpfig(gdata(gname),ncpdata(gname).ncp,0,0,0);
    print(sprintf('nvsncp-%s.eps',gname),'-depsc2','-painters');
end
