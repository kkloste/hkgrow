% run from /kdd
% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convtwit2symmat > convtwit.txt &

filename = 'twitter-2010';
P = load_graph(filename,'/scratch/dgleich/kyle/data');
n = size(P,1);
P = P|P;

outputname = strcat('twitter');
save(['/scratch/dgleich/kyle/symmats/' outputname '.mat'], 'P', 'n','-v7.3');
exit;

