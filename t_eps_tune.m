% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r t_eps_tune > /dev/null 2>&1&
% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r t_eps_tune > tepstune.txt &


load /scratch2/dgleich/kyle/symmats/ljournal;
filename = 'ljournal';

% setup inputs
% A = load_graph(filename,'~/data'); n = size(A,1);


% numeps = number of different values for epsilon
% numt = number of different values for t
% numtrials = number of trials

debugflag = 0;

numtrials = 10;

n = size(P,1);

eps_vals = [ 1e-2 5*1e-3 1e-3 5*1e-4 1e-4]';
t_vals = [30 70 80 100]';
numt = numel(t_vals);
numeps = numel(eps_vals);


targetvol = 1000;
indices = randi(n,numtrials,1);
times = zeros(numeps,numt,numtrials);
conds = zeros(numeps,numt,numtrials);

% Test on random seeds
for trial_num=1:numtrials
fprintf( 'trial= %i \t', trial_num);
for eps_num=1:numeps
    for t_num=1:numt
tic; [dummy conds(eps_num,t_num,trial_num) cut vol] = hkgrow_mex(P,indices(trial_num),targetvol,t_vals(t_num), eps_vals(eps_num), debugflag);
% tic; [dummy, conds(eps_num,t_num,trial_num), cut, vol] = hkgrow(P,indices(trial_num),eps_vals(eps_num), t_vals(t_num), debugflag);

        times(eps_num,t_num,trial_num) = toc;
    end
end
end

avconds = zeros(numeps,numt);
avtimes = avconds;
for i = 1:numeps, for j=1:numt
    avconds(i,j) = sum(conds(i,j,:));
    avtimes(i,j) = sum(times(i,j,:));
end, end
avconds = avconds./numtrials
avtimes = avtimes./numtrials

outputname = strcat(filename,'GRIDteps');
save(['/scratch2/dgleich/kyle/results/' outputname '.mat'], 'avtimes', 'avconds', 'eps_vals', 't_vals' ,'-v7.3');
% matlabmail hasn't been working, comment out for now
% matlabmail('kyle.kloster@gmail.com', 'experiment [t_eps_tune] done', 'matlab - [t_eps_tune] done');
exit;
