% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r experimenthkgrow > exphk.txt &

filename={'pgp-cc', 'ca-AstroPh-cc', 'marvel-comics-cc', 'as-22july06', 'rand-ff-25000-0.4', 'cond-mat-2003-cc', 'email-Enron-cc','cond-mat-2005-fix-cc', 'soc-sign-epinions-cc', 'itdk0304-cc', 'dblp-cc', 'flickr-bidir-cc'};
numtrials = 200;

indices = zeros(numel(filename),numtrials,4);
times = zeros(numel(filename),numtrials,4);
conds = zeros(numel(filename),numtrials,4);
gsize = zeros(numel(filename),2);

for fileid=1:numel(filename)
    dataset = char(filename(fileid));
    fprintf('file number = %i  ,  graph = %s \n', fileid, dataset);
    A = load_graph(dataset,'/scratch2/dgleich/kyle/data');
%    A = A|A';
    gsize(fileid,1) = size(A,1);
    gsize(fileid,2) = nnz(A);

    [conds(fileid,:,1) times(fileid,:,1) indices(fileid,:,1)] = randseed(A,numtrials);
fprintf('\t randseed done\n');
    [conds(fileid,:,2) times(fileid,:,2) indices(fileid,:,2)] = heavyseed(A,numtrials);
fprintf('\t heavyseed done\n');
    [conds(fileid,:,3) times(fileid,:,3) indices(fileid,:,3)] = randhood(A,numtrials);
fprintf('\t randhood done\n');
    [conds(fileid,:,4) times(fileid,:,4) indices(fileid,:,4)] = heavyhood(A,numtrials);
fprintf('\t heavyhood done\n');
    clear A;
end

outputname = strcat('smalldata');
save(['/scratch2/dgleich/kyle/results/' outputname '.mat'], 'gsize', 'indices', 'times', 'conds', 'filename','-v7.3');
clear
exit;


