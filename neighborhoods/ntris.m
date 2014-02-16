function t = ntris(A)
% NTRIS Compute the number of triangles centered at each vertex of a graph.
%
% This is just a convinence wrapper around
%  d = sum(A,2);
%  t = clustering_coefficients(A)*(d-1)*(d)/2;
%  

assert(isequal(A,A'));

d = full(sum(A,2));
t = clustering_coefficients(A);
t = t.*d.*(d-1)/2;
t = round(t); % make sure we get out integers
