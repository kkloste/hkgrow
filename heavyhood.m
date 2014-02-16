function [conds times indices] = heavyhood(A,numtrials,tol,alphat)
% [conds times indices] = heavyhood(A,numtrials,tol,alphat)

n = size(A,1);
degrees = zeros(n,1);
for ind = 1:n
    degrees(ind) = nnz(A(:,ind));
end
[vals, perm] = sort(degrees,'descend');
indices = perm(1:numtrials);
 
times = zeros(numtrials,1);
conds = zeros(numtrials,1);
 
for trial_num=1:numtrials
    [neighborhood, ~, ~] = find(A(:,indices(trial_num)));
    tic; [dummy,conds(trial_num),cut_hk,vol_hk] = hkgrow(A,neighborhood,tol,alphat,0);
    times(trial_num) = toc;
end
