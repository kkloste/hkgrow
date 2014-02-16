function cond = cutcond(A,s)
% CUTCOND Return the conductance of a cut
%
% cond = cutcond(A,s) returns the sum of degrees of vertices in A
%

d = sum(A,2);
Gvol = full(sum(d));
setvol = cutvol(A,s);
cut = cutsize(A,s);
cond = cut./min(Gvol-setvol,setvol);