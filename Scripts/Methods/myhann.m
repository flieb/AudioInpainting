function g = myhann(N,T,L)

g = zeros(L,1);
hann = @(x) 0.5*(1-cos(2*pi*(0:x-1)/(x-1)));
tmp = hann(N);
g(1+T:T+N-1) = tmp(1:end-1);