%%
eps = 5*1e-3;
t = 1;
k=0;
while k < 250,
    epsexpt = exp(t)*eps;
    error = exp(t) - 1;
    last = 1;
    k = 0;
    while error > epsexpt,
        k = k+1;
        last = (last*t)/k;
        error = error - last;
    end
   [k t] 
   t=t+1;
end