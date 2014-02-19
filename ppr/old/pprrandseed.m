function [conds times indices] = pprrandseed(A,numtrials)
% [conds times indices] = pprrandseed(A,numtrials)

n = size(A,1);

indices = randi(n,numtrials,1);
times = zeros(numtrials,1);
conds = zeros(numtrials,1);

for trial_num=1:numtrials
    tic; [dummy,conds(trial_num),cut_hk,vol_hk] = pprgrow(A,indices(trial_num));
    times(trial_num) = toc;
end
