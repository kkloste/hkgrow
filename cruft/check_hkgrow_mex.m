 A = load_graph('netscience-cc','~/data');
% A = load_graph('email','~/data'); %just trying another, to be sure

n = size(A,1);
c = randi(n);
t = 5;
eps = 1e-3;
debugflag = 0;

tic; [bset, cond, cut, vol, hkvec, npush] = hkgrow_mex(A,c,t,eps,debugflag);
toc
