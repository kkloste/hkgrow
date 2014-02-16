function [cond cut vol csize] = ball_clusters(A,k)
n = size(A,1);
Ak = A;
Ag = A;
for i=2:k
    Ak = Ak*A;
    Ag = Ak + Ag;
end
Ag = spones(Ag + speye(n));

cond = zeros(n,1);
cut = zeros(n,1);
vol = zeros(n,1);
csize = zeros(n,1);

for i=1:n
    set = find(Ag(:,i));
    [cond(i) cut(i) vol(i)] = cutcond_mex(A,set);
    csize(i) = min(n-numel(set),numel(set));
end