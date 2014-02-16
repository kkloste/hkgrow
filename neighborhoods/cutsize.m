function cs=cutsize(A,s)
% Only works on symmetric graphs

n = size(A,1);

svec = zeros(n,1);
svec(s) = 1;

if length(svec)>n/2
    svec = 1-svec;
    s = find(svec);
end

[i,j] = find(A(:,s));
cs = 0;
nzs = length(i);
for nzi=1:nzs
    if svec(i(nzi))==0
        cs = cs + 1;
    end
end
