% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r trials_twitter > trialstwitter.txt &

load /scratch2/dgleich/kyle/symmats/twitter;
filename = 'twitter';

numtrials = 10;
tvals = 55;
eps = 1e-3;
indices = zeros(numtrials,1);
times = zeros(numtrials,1);
conds = zeros(numtrials,1);
gsize = zeros(1,2);

fprintf('starting graph = %s\n', filename);
gsize(1) = size(A,1);
gsize(2) = nnz(A);
[conds(:,1) times(:,1) indices(:,1)] = randseed(A,numtrials,eps,tvals);
fprintf('randseed done  \n');
%    [conds(:,1) times(:,1) indices(:,1)] = heavyseed(A,numtrials,eps,tvals);
%    fprintf('heavyseed done  \n');
%    [conds(:,1) times(:,1) indices(:,1)] = randhood(A,numtrials,eps,tvals);
%    fprintf('randhood done  \n');
%    [conds(:,1) times(:,1) indices(:,1)] = heavyhood(A,numtrials,eps,tvals);

avecond = sum(conds(:,1))./numtrials;
avetime = sum(times(:,1))./numtrials;
fprintf('avecond=%f  avetime=%f \n', avecond, avetime);
outputname = strcat(filename,'trials');
save(['/scratch2/dgleich/kyle/results/' outputname '.mat'], 'tvals','eps','gsize', 'indices', 'times', 'conds', 'filename','-v7.3');
% matlabmail('kyle.kloster@gmail.com', 'experiment [trials_twitter] done', 'matlab - [trials_twitter] done');
exit;