BCtest = sparse([0 1 1 0 0; 1 0 1 0 0; 1 1 0 1 0; 0 0 1 0 1; 0 0 0 1 0]);

[BCcc BCp BCf BCfcc] = biconncore(BCtest);

assert(all(all(BCcc == sparse([0 1 1; 1 0 1; 1 1 0]))));
assert(all(BCp == [1 1 1 0 0]'));
assert(all(all(BCf == ...
    sparse([0 1 1 0 0; 1 0 1 0 0; 1 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0]))));
assert(max(BCfcc) == 2);