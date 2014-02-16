%% Plots of the properties of the neighborhood communities.
load ../../data/gdata

%% Show them all
% Here, we'll pick out a few interesting figures to put into
% the paper based on properties of all of them.

for gname=gdata.keys
    g = gname{1};
    [cond i] = min(gdata(g).neigh.cond);
    fset = gdata(g).fiedler.set;
    i_in_f = any(fset==i);
    fprintf('%35s : %i in f   %.3f mincond  %.3f fcond\n',g,...
        i_in_f,min(gdata(g).neigh.cond(fset)),gdata(g).fiedler.cond);
end 

