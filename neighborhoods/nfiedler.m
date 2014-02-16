function [x,p,cuts,bestset,lam2] = nfiedler(A,tol)
% NFIEDLER Return the normalized Fiedler vector for a graph
%
% [x,p,cuts,bestset,lam2] = nfielder(A);

if nargin<2, tol=1e-5; end

L = nlaplacian(A);
n = size(A,1);
[V,D] = eigs(L+speye(n),2,'SA',struct('tol',tol));
x = V(:,2);
x = x./sqrt(sum(A,2));
lam2 = D(2,2)-1;
if nargout>1
    [~,p] = sort(x);
end
if nargout>2
    cuts = cutsweep(A,p);
    [val,ind] = min(cuts.conductance);
    if ind<((n+1)/2)
        bestset = cuts.order(1:ind);
    else
        bestset = cuts.order(ind+1:end);
    end
end

