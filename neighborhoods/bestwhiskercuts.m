function [cond,cut,vol,s,p] = bestwhiskercuts(A)

[~,~,volw,sw,fw] = whiskercuts(A);
[~,pw] = sort(volw,'descend');
Gvol = sum(A(:));
nwhis = length(pw);
cond = zeros(nwhis,1);
cut = zeros(nwhis,1);
vol = zeros(nwhis,1);
s = zeros(nwhis,1);
p = Inf*ones(size(A,1),1);

cursize = 0;
curvol = 0;
for i=1:length(pw)
    wi = pw(i);
    p(fw==wi) = i; % these are in the ith cluster
    
    cursize = cursize + sw(wi);
    curvol = curvol + volw(wi);
    
    cut(i) = i;
    s(i) = cursize;
    vol(i) = curvol;
    cond(i) = i./min(curvol,Gvol-curvol);
end
    
    
    