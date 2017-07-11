function [g,BL,BR,l] = getwin(desired_support, bins)


% gets wp2inp of correct support size

if desired_support < 3
    error('support is lower than 3, not possible... check code');
end

a = 2^(1/bins);

supp = 0;

load linefit.mat

j = floor(desired_support*linefit(1) + linefit(2));
ds = floor(desired_support);

while supp ~= ds

    dx = linspace(1e-3,3,j);
    x = log(dx)/log(a);
    mu = 25;

    g = exp(a.^(x).*(1-log(a.^(x)))).^mu;
    %g = exp(exp(-2*x)*mu.*(1+2*x));
    
    %tic,
    [BL,BR,supp] = support(g,1e-7);
    %toc
    
    
    %fprintf('supp = %d\n',supp);
    
    if supp < ds
        j = j+8;
    elseif supp > ds
        j = j-1;
    end
   
end
g = g'./max(g);
l = j;