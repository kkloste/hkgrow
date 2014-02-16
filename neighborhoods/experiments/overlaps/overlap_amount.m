%% Check the number of good clusters we can accumulate without overlap

%%
A = load_external_graph('snap/email-Enron');
%A = load_external_graph('snap/cit-HepTh');
%A = load_graph('as-22july06');
%A = load_graph('cond-mat-2005-fix-cc');
n = size(A,1);

%%
cond = triangleclusters_mex(A);
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
[~,p] = sort(cond./d);

nassigned=0;
assigned=false(n,1);
picked=false(n,1);

sumcond = 0;

for i=1:n
    v = p(i);
    nv = find(A(:,v));
    if assigned(v) || any(assigned(nv))
        continue
    end
    
    if cond(v) > worstcond
        continue;
    end
    
    assigned(v) = true;
    assigned(nv) = true;
    nassigned = nassigned + length(nv);
    if ~A(v,v), nassigned = nassigned+1; end
    
    picked(v) = true;
    
    sumcond = sumcond + cond(p(i));
end

clf;
ncpplot(cond,d+1); hold on;
plot(d(picked)+1,cond(picked),'ro');
hold off;
title(sprintf('%i/%i assigned   %f mean cond   %f max cond', ...
    sum(assigned), n, mean(cond(picked)), worstcond));
%%
% What's left should be really dull
Iden = speye(n);
AI = A-diag(diag(A))+Iden;
R = [AI(:,picked) Iden(:,~assigned)];
assert(n-sum(assigned)+sum(picked) == size(R,2));
B = spones(R'*A*R);
B = B - diag(diag(B));
[f,fp,cuts,bestset,lam2] = nfiedler(B);