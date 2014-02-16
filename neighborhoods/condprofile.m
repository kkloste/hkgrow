function h1 = condprofile(cond,size)
A = accumarray(size,cond,[], @condhist, {{[],[]}});

% extract the data
npts = sum(cellfun(@(c) numel(c{1}),A));
x = zeros(npts,1);
y = zeros(npts,1);
cnt = zeros(npts,1);nz = cnt>0;
scatter(x(nz),y(nz),log10(cnt(nz)));
curpt = 1;
for i=1:numel(A)
    nbins = numel(A{i}{1});
    x(curpt:curpt+nbins-1) = i*ones(nbins,1);
    y(curpt:curpt+nbins-1) = 10.^(A{i}{2}); % the points x
    cnt(curpt:curpt+nbins-1) = A{i}{1};
    curpt = curpt + nbins;
end
assert(curpt == npts+1);
nz = cnt>0;
nclus = sum(cnt(nz));
scatter(x(nz),y(nz),10000*cnt(nz)./nclus);
set(gca,'XScale','log');
set(gca,'YScale','log');

function celldat = condhist(cond)
if isempty(cond), celldat = {[],[]}; return; end
nbins = floor(5*(-log10(min(cond))));
[n,x] = hist(log10(cond),nbins);
celldat = {{n,x}};

