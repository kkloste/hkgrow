function [conds times indices] = randhood(A,numtrials,tol,alphat)
% [conds times indices] = randhood(A,numtrials,tol,alphat)

n = size(A,1);

indices = randi(n,numtrials,1);
times = zeros(numtrials,1);
conds = zeros(numtrials,1);

for trial_num=1:numtrials
    [neighborhood, ~, ~] = find(A(:,indices(trial_num)));
    tic; [dummy,conds(trial_num),cut_hk,vol_hk] = hkgrow(A,neighborhood,tol,alphat,0);
    times(trial_num) = toc;
end
