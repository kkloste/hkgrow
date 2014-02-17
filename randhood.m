function [conds times indices] = randhood(A,numtrials)
% [conds times indices] = randhood(A,numtrials)

n = size(A,1);

indices = randi(n,numtrials,1);
times = zeros(numtrials,1);
conds = zeros(numtrials,1);

for trial_num=1:numtrials
    tic; [dummy,conds(trial_num),cut_hk,vol_hk] = hkgrow(A,indices(trial_num),'neighborhood',true);
    times(trial_num) = toc;
end
