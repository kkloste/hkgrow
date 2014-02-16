 A = load_graph('netscience-cc','~/data');
% A = load_graph('email','~/data'); %just trying another, to be sure
ntrials = 50;
for ni = 1:size(A,1);
    for tols = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7]
        [set1,cond1] = hkgrow_sresid_mex(A,ni,100000,1,tols,0);
        [set2,cond2] = hkgrow_mex(A,ni,100000,1,tols,0);
        assert(cond1 == cond2);
    end
end