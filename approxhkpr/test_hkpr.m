function [time_hk cond_hk bestset_hk setup_time] = test_hkpr(filename,numtrials,tol,alphat,debugflag)
% [times conductances cut_sets setup_time] = test_hkpr(filename,numtrials,tol,alphat,debugflag)
% set debugflag to 1 to turn on messages in the matlab and mex code.
tic;

if nargin < 1,
filename = 'netscience-cc';
end

if nargin < 2,
numtrials = 1;
end

if nargin < 3,
tol  = 1e-5;
end

if nargin < 4,
alphat = 1;
end

if nargin < 5,
debugflag = 0;
end

assert(tol > 0 && tol <= 1, 'tol violates 0<tol<=1');
assert(alphat>0, 'alphat violates alphat>0');
assert(numtrials>=1, 'numtrials must be positive integer');


A = load_graph(filename,'~/data'); n = size(A,1);
degrees = zeros(n,1);
for ind = 1:n
degrees(ind) = nnz(A(:,ind));
end
setup_time = toc;

if debugflag == 1,
fprintf('test_hkpr: setup time=%f \n', setup_time);
end

%% First do random seeds
indices = randi(n,numtrials,1);

time_hk = zeros(numtrials,1);
cond_hk = zeros(numtrials,1);
bestset_hk = zeros(n,numtrials);

for trial_num=1:numtrials

if debugflag==1, fprintf('test_hkpr:  start rand trial=%i \n', trial_num); end

tic; [dummy,cond_hk(trial_num),cut_hk,vol_hk] = hkprgrow(A,indices(trial_num),tol,alphat,debugflag);

if debugflag == 1, fprintf('test_hkpr:  end rand trial=%i \n', trial_num); end

bestset_hk(1:size(dummy,1),trial_num) = dummy;
time_hk(trial_num) = toc;
end

