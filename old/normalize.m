function [A,d,dmax] = normalize(A,n)
% [A,d] = NORAMLIZE(A,n)
% updates A in place (to reduce memory use)
% with a column stochastic version of A

d = zeros(n,1); dmax = 0;
for ind=1:n
    d(ind) = nnz(A(:,ind));
    A(:,ind) = A(:,ind)./d(ind);
    if (dmax < d(ind))
        dmax = d(ind);
    end
end
