function vol = cutvol(A,s)
% CUTVOL Return the volume of a cut
%
% vol = cutvol(A,s) returns the sum of degrees of vertices in A
%

d = sum(A,2);
vol = full(sum(d(s)));