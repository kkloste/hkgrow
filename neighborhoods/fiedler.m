function [x,p] = fiedler(A)
% FIEDLER Return the Fiedler vector for a graph

L = laplacian(A);
[V,D] = eigs(L,2,'SA');
x = V(:,2);
[~,p] = sort(x);
