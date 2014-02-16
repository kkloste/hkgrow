addpath('../../matlab/');
A = readSMAT('../../data/polblogs-sym-cc.smat');
H = multicluster(A,[1],'maxcond',1);
 