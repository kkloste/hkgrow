

function best_cond = trials_hkgrow(filename,numtrials,saveasname)
% best_cond = trials_pprgrow(filename,numtrials,saveasname)

% load /scratch/dgleich/kyle/mat-files/[name of dataset];
cd ~/kdd;
A = load_graph(filename,'~/data'); P = colnormout(A); clear A; n = size(P,1);
cd ~/kdd/neighborhoods;

if (nargin<2)
numtrials = 100;
end

%% First do random seeds

% load the indices from the hkgrow trials
load(strcat('~/kdd/',filename,'-',saveasname,'-randseed-hkgrow'));

time_ppr = zeros(numtrials,1);
cond_ppr = zeros(numtrials,1);
bestset_ppr = zeros(n,numtrials);
for trial_num=1:numtrials
tic; [dummy,cond_ppr(trial_num),cut_ppr,vol_ppr] = pprgrow(P,indices(trial_num)); time_ppr(trial_num) = toc;
bestset_ppr(1:size(dummy,1),trial_num) = dummy;

end

save(strcat(filename,'-',saveasname,'-randseed-pprgrow'), 'indices', 'time_ppr', 'cond_ppr', 'bestset_ppr' ,'-v7.3');



%% Now do degree sorted seeds

% load the indices from the hkgrow trials
load(strcat('~/kdd/',filename,'-',saveasname,'-heavyseed-hkgrow'));

time_ppr = zeros(numtrials,1);
cond_ppr = zeros(numtrials,1);
bestset_ppr = zeros(n,numtrials);
for trial_num=1:numtrials
tic; [dummy,cond_ppr(trial_num),cut_ppr,vol_ppr] = pprgrow(P,indices(trial_num)); time_ppr(trial_num) = toc;
bestset_ppr(1:size(dummy,1),trial_num) = dummy;
time_hk(trial_num) = toc;
end

save(strcat(filename,'-',saveasname,'-heavyseed-pprgrow'), 'indices', 'time_ppr', 'cond_ppr', 'bestset_ppr' ,'-v7.3');
