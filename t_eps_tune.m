% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r t_eps_tune > /dev/null 2>&1&
% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r t_eps_tune > results.txt &

tic;
load /scratch/dgleich/kyle/mat-files/ljournal-2008;
filename = 'ljournal';
A = spones(P);
A = A|A;
clear P;
n = size(A,1);

%assert(filename=!NULL, 'enter a filename');

% setup inputs
% A = load_graph(filename,'~/data'); n = size(A,1);
setup_time = toc;

% numeps = number of different values for epsilon
% numt = number of different values for t
% numtrials = number of trials

debugflag = 0;

numtrials = 100;

% eps_vals = zeros(numeps,1);
% t_vals = zeros(numt,1);

eps_vals = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6]';
t_vals = [1 5 15 20 25 30]';
numt = numel(t_vals);
numeps = numel(eps_vals);


curexpand = 1000;
indices = randi(n,numtrials,1);
time_hk = zeros(numeps,numt,numtrials);
cond_hk = zeros(numeps,numt,numtrials);
vol_hk = zeros(numeps,numt,numtrials);

% Test on random seeds
for trial_num=1:numtrials
fprintf( 'trial= %i  : eps = %i  t = %i\n', trial_num, eps_num, t_num);
for eps_num=1:numeps
    for t_num=1:numt
fprintf( 'eps = %i  t = %i \t', trial_num, eps_num, t_num);
        tic; [dummy cond_hk(eps_num,t_num,trial_num) cut_hk vol_hk(eps_num,t_num,trial_num)] = hkgrow_mex(A,indices(trial_num),curexpand,t_vals(t_num), eps_vals(eps_num), debugflag);
        time_hk(eps_num,t_num,trial_num) = toc;
% kappai version
%tic; [dummy cond_hk(eps_num,t_num,trial_num,2) cut_hk vol_hk(eps_num,t_num,trial_num,2)] = kap_hkgrow_mex(A,indices(trial_num),curexpand,% t_vals(t_num), eps_vals(eps_num), debugflag);
%        time_hk(eps_num,t_num,trial_num,2) = toc;
    end
end
end

outputname = strcat(filename,'GRIDtwitter');
save(['/scratch/dgleich/kyle/results/' outputname '.mat'], 'indices', 'time_hk', 'cond_hk', 'eps_vals', 't_vals' ,'-v7.3');

exit;
