function [n] = myl21norm(x, param)

%computes the l21 norm 

n = norm(computeWeight(x,param),1);