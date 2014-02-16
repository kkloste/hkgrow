% run from /kdd
% nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convtwit2symmat > convtwit.txt &

filename = 'com-friendster.ungraph';
P = load_graph(filename,'/scratch/dgleich/friendster');
n = size(P,1);
    fprintf('done loading, now symmetrizing \n');
P = P|P;

outputname = strcat('friendster');
save(['/scratch/dgleich/kyle/symmats/' outputname '.mat'], 'P', 'n','-v7.3');
exit;