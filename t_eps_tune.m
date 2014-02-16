% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r t_eps_tune > /dev/null 2>&1&
% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r t_eps_tune > jourgrid.txt &


load /scratch/dgleich/kyle/symmats/ljournal;
filename = 'ljournal';

% setup inputs
% A = load_graph(filename,'~/data'); n = size(A,1);


% numeps = number of different values for epsilon
% numt = number of different values for t
% numtrials = number of trials

debugflag = 0;

numtrials = 100;

eps_vals = [5*1e-2 1e-2 1e-3 1e-4 1e-5]';
t_vals = [1 5 15 20 25 30]';
numt = numel(t_vals);
numeps = numel(eps_vals);


curexpand = 1000;
indices = randi(n,numtrials,1);
times = zeros(numeps,numt,numtrials);
conds = zeros(numeps,numt,numtrials);

% Test on random seeds
for trial_num=1:numtrials
fprintf( 'trial= %i \t', trial_num);
for eps_num=1:numeps
    for t_num=1:numt
        tic; [dummy conds(eps_num,t_num,trial_num) cut vol] = hkgrow_sresid_mex(P,indices(trial_num),curexpand,t_vals(t_num), eps_vals(eps_num), debugflag);
        times(eps_num,t_num,trial_num) = toc;
    end
end
end

outputname = strcat(filename,'GRIDteps');
save(['/scratch/dgleich/kyle/results/' outputname '.mat'], 'indices', 'times', 'conds', 'eps_vals', 't_vals' ,'-v7.3');

exit;
