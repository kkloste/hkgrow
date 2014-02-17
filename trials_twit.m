% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r trials_twit > trialstwit.txt &

load /scratch/dgleich/kyle/symmats/twitter;
filename = 'twitter';

numtrials = 50;
tvals = 20;
eps = 5*1e-3;
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