function [fun,k1,kend,j] = mwavelet(name,a,fmin,bw,fs,epsilon)

%returns the anonymous wavelet supported on the interval [M1,M2] such that
%1-eps are containted in this interval. 

if nargin < 5
    epsilon = 1e-6;
end

W = @(x,k,j) j*log(x)./log(a) - j*k;

k1 = log(fmin)./log(a);
b = bw/(a^k1);

switch name
    case 'gaussian'
        ogauss = @(x,sigma) exp(-( (a.^(x)-5).^2)/(2*sigma));
        
        
        k1 = log(fmin)/log(a) - log(5)/log(a);
        sigma = (a^(-k1)*bw/2/erfinv(1-epsilon))^2;
        %M2 = a^k1*sqrt(sigma)*erfinv((1-epsilon)) + a^k1*5;
        j = 1;
        kend = log((fs/2)/(erfinv(1-epsilon)*sqrt(sigma) + 5))/log(a);
        
%         x = -7:0.0001:7;
%         sigma = 2;
%         g = ogauss(x,sigma);
%         [l,r] = support(g,epsilon);
%         M1 = x(l);
%         M2 = x(r);
%         fun = @(x) x.^M2 - x.^M1 - b;
%         y = fzero(fun,[epsilon 100]);
%         j = log(a)/log(y);
%         kend = log(fs/2)/log(a) - M2/j;
        
        
        
        fun = @(x,k,j) ogauss(W(x,k,j),sigma);
    case 'warpedGaussian'
        gauss = @(x) exp(-( (x).^2)/2);
        x = -7:0.0001:7;
        g = gauss(x);
        [l,r] = support(g,epsilon);
        M1 = x(l);
        M2 = x(r);
        
        fun = @(x) x.^M2 - x.^M1 - b;
        y = fzero(fun,[epsilon 100]);
        j = log(a)/log(y);
        
        kend = log(fs/2)/log(a) - M2/j;
        
        fun = @(x,k,j) exp(-(W(x,k,j).^2)/2);
        
    case 'wp2inp'
        wpinp = @(x) 1/exp(25)*(a.^x>=0).*exp(25*a.^x.*(1-log(a.^x)));
        x = -2:0.00001:2;
        g = wpinp(x);
        [l,r] = support(g,epsilon);
        M1 = x(l); 
        M2 = x(r);
        
        fun = @(x) x.^M2 - x.^M1 - b;
        y = fzero(fun,[epsilon 100]);
        j = log(a)/log(y);
        
        kend = log(fs/2)/log(a) - M2/j;
        
        fun = @(x,k,j) 1/exp(25)*(a.^(W(x,k,j))>=0).*exp(25*(a.^W(x,k,j)).*(1-log(a.^W(x,k,j))));
    case 'warpedWp2inp'
        wpinp = @(x) 1/exp(25)*(x+1>=0).*exp(25*(x+1).*(1-log(x+1)));
        x = (-1 + 0.00001):0.00001:10;
        g = wpinp(x);
        [l,r] = support(g,epsilon);
        M1 = x(l);
        M2 = x(r);
        
        fun = @(x) x.^M2 - x.^M1 - b;
        y = fzero(fun,[epsilon 100]);
        j = log(a)/log(y);
        
        kend = log(fs/2)/log(a) - M2/j;
        
        fun = @(x,k,j) 1/exp(25)*((W(x,k,j))+1>=0).*exp(25*(W(x,k,j)+1).*(1-log(W(x,k,j)+1)));
        
    case 'mexicanHat'
        fmxh = @(x) 8*sqrt(2/3)*pi^(9/4)*(a.^(x.^2)).*exp(-2*pi^2*a.^(x.^2));
        x = -1:0.001:4;
        g = fmxh(x);
        
        [l,r] = support(g,epsilon);
        M1 = x(l);
        M2 = x(r);
        
        fun = @(x) x.^M2 - x.^M1 - b;
        y = fzero(fun,[epsilon 100]);
        j = log(a)/log(y);
        
        kend = log(fs/2)/log(a) - M2/j;
        fun = @(x,k,j) 8*sqrt(2/3)*pi^(9/4)*(a.^(W(x,k,j))>=0).*(a.^(W(x,k,j).^2)).*exp(-2*pi^2*a.^(W(x,k,j).^2));
    otherwise
        error('not defined');
end

