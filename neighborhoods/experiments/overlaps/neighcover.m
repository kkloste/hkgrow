function [maxconds,coverage,sumconds] = neighcover(A, cond,p)
% NEIGHCOVER Cover the grah with neighborhoods
%
% cond = conductance of each neighborhood
% p order in which to examine neighborhoods

n = size(A,1);

nassigned=0;
assigned=false(n,1);
picked=false(n,1);
npicked=0;

sumcond = 0;
maxcond = -1;

maxconds = zeros(n,1);
coverage = zeros(n,1);
sumconds = zeros(n,1);

for i=1:n
    v = p(i);
    nv = find(A(:,v));
    if assigned(v) || any(assigned(nv))
        continue
    end
    
    assigned(v) = true;
    assigned(nv) = true;
    nassigned = nassigned + length(nv)+1;
    
    
    picked(v) = true;
    npicked = npicked + 1;
    
    sumcond = sumcond + cond(p(i));
    maxcond = max(cond(p(i)),maxcond);
    maxconds(npicked) = maxcond;
    coverage(npicked) = nassigned/n;
    sumconds(npicked) = sumcond;
end

maxconds = maxconds(1:npicked);
coverage = coverage(1:npicked);
sumconds = sumconds(1:npicked);