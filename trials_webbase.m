% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r experimenthkgrow > hkgrowsmallsets.txt &

load /scratch/dgleich/kyle/symmats/webbase;
filename = 'webbase';

numtrials = 100;
tvals = 20;
eps = 1e-3;
indices = zeros(numel(filename),numtrials,4);
times = zeros(numel(filename),numtrials,4);
conds = zeros(numel(filename),numtrials,4);
gsize = zeros(numel(filename),2);
for fileid=1:numel(filename)
    dataset = char(filename(fileid));
    fprintf('file number = %i  ,  graph = %s\n', fileid, dataset);
gsize(fileid,1) = size(P,1);
gsize(fileid,2) = nnz(P);
    [conds(fileid,:,1) times(fileid,:,1) indices(fileid,:,1)] = randseed(P,numtrials,eps,tvals);
    fprintf('randseed done  \n');
    [conds(fileid,:,2) times(fileid,:,2) indices(fileid,:,2)] = heavyseed(P,numtrials,eps,tvals);
    fprintf('heavyseed done  \n');
    [conds(fileid,:,3) times(fileid,:,3) indices(fileid,:,3)] = randhood(P,numtrials,eps,tvals);
    fprintf('randhood done  \n');
    [conds(fileid,:,4) times(fileid,:,4) indices(fileid,:,4)] = heavyhood(P,numtrials,eps,tvals);
end
datalabels = 'indices(fileid,numtrial,experimenttype)';
outputname = strcat(filename,'hkgrow');
save(['/scratch/dgleich/kyle/results/' outputname '.mat'], 'tvals','eps','gsize', 'indices', 'times', 'conds', 'filename','-v7.3');
exit;






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




