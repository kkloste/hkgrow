function [conds times indices] = randseed(A,numtrials,tol,alphat)
% [conds times indices] = randseed(A,numtrials,tol,alphat)

n = size(A,1);

indices = randi(n,numtrials,1);
times = zeros(numtrials,1);
conds = zeros(numtrials,1);

for trial_num=1:numtrials
    tic; [dummy,conds(trial_num),cut_hk,vol_hk] = hkgrow(A,indices(trial_num),tol,alphat,0);
    times(trial_num) = toc;
end
