%% Check the number of good clusters we can accumulate without overlap

%%
A = load_external_graph('snap/email-Enron');
%A = load_external_graph('snap/cit-HepTh');
%A = load_graph('as-22july06');
%A = load_graph('cond-mat-2005-fix-cc');
n = size(A,1);

%%
[cond,cut,vol,s,t,cc,Cdata] = triangleclustersgreedy_mex(A);
C = sparse(Cdata(:,2),Cdata(:,1),1,n,n);
d = full(sum(A,2));
%%
%cdplot(cond,d+1);
%%
bestcond = min(cond);
%worstcond = mean(cond);
%worstcond = mean(cond);
worstcond = prctile(cond(d>1),10);
%worstcond = mean(cond);
% sort by cond
[~,p] = sort(cond./s);

nassigned=0;
assigned=false(n,1);
picked=false(n,1);

sumcond = 0;

for i=1:n
    v = p(i);
    nv = find(C(:,v));
    if assigned(v) || any(assigned(nv))
        continue
    end
    
    if cond(v) > worstcond
        continue;
    end
    
    assigned(v) = true;
    assigned(nv) = true;
    nassigned = nassigned + length(nv);
    
    picked(v) = true;
    
    sumcond = sumcond + cond(p(i));
end

clf;
ncpplot(cond,s+1); hold on;
plot(s(picked)+1,cond(picked),'ro');
hold off;
title(sprintf('%i/%i assigned   %f mean cond   %f max cond', ...
    sum(assigned), n, mean(cond(picked)), worstcond));
