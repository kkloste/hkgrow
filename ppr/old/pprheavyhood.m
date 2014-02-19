function [conds times indices] = pprheavyhood(A,numtrials)
% [conds times indices] = pprheavyhood(A,numtrials)

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
    tic; [dummy,conds(trial_num),cut_hk,vol_hk] = pprgrow(A,indices(trial_num),'neighborhood',true);
    times(trial_num) = toc;
end
