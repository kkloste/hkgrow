%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
A = load_graph('email','~/data'); n = size(A,1);
P = normout(A); clear A;
num_test = 5;

x_true = zeros(n,num_test);
time_vals = zeros(num_test,2);
error_vals = zeros(num_test,2);
set_vals = zeros(num_test,2,2);

%
% times_vals[i,j] = trial i, method j
% j : 1 = true , 2 = approx_hkpr
%

tol = 1e-1;
k=[10 100];

for i=1:num_test
i = i
eyei = zeros(n,1);
c = randi(n);
eyei(c)=1;

tic; [x_true(:,i),sdummy,m,mv,mvd] = expmv(1,P,eyei,[],'single'); time_vals(i,1) = toc; normtrue = norm(x_true(:,i),1); [vtrue strue] = sort(x_true(:,i), 'descend');

tic; [x_h n_pushes] = approxhkpr_mex(P,c,1,tol,debugflag); time_vals(i,2) = toc;

error_vals(i,2) = norm(x_h - x_true(:,i),1)/normtrue;


[vh sh1] = sort(x_h, 'descend');
%[vh sh1] = sort(x_gexpmq, 'descend');
%[vh sh2] = sort(x_gexpm, 'descend');


for t = 1:2
set_vals(i,2,t) = length(intersect(sh1(1:k(t)),strue(1:k(t))))/k(t);
end
end



