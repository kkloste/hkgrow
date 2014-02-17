% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r trials_ljournal > trialsljournal.txt &

load /scratch2/dgleich/kyle/symmats/ljournal;
filename = 'ljournal';

numtrials = 100;
tvals = 25;
eps = 1e-3;
indices = zeros(numtrials,2);
times = zeros(numtrials,2);
conds = zeros(numtrials,2);
gsize = zeros(1,2);

    fprintf('starting graph = %s\n', filename);
gsize(1) = size(P,1);
gsize(2) = nnz(P);
    [conds(:,1) times(:,1) indices(:,1)] = randseed(P,numtrials,eps,tvals);
    fprintf('randseed done  \n');
    [conds(:,2) times(:,2) indices(:,2)] = heavyseed(P,numtrials,eps,tvals);
    fprintf('heavyseed done  \n');
%    [conds(:,1) times(:,1) indices(:,1)] = randhood(P,numtrials,eps,tvals);
%    fprintf('randhood done  \n');
%    [conds(:,1) times(:,1) indices(:,1)] = heavyhood(P,numtrials,eps,tvals);

outputname = strcat(filename,'trials');
save(['/scratch2/dgleich/kyle/results/' outputname '.mat'], 'tvals','eps','gsize', 'indices', 'times', 'conds', 'filename','-v7.3');
% matlabmail didn't work here, can't explain why, comment out for now
% matlabmail('kyle.kloster@gmail.com', 'experiment [trials_ljournal] done', 'matlab - [trials_ljournal] done');
exit;