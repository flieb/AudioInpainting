function [sol,W] = myerbproxl21(x,lambda,param)

%l21 proximal operator for nonstationary transforms
%input signal is a concatenated vector of all data points

M = param.M;
l = param.t_length;
L = param.L;
fs = param.fs;
cM = [1; M];
cM = cumsum(cM);

st = @(s,T,W) s.*max(1-T./abs(W),0);
pew = @(s,T,W) s.*max(1-(T./abs(W)).^2,0);

W = x;


nkernel = max(1,round( l./(L*1000./(fs.*M) )));
indx = find(nkernel-1);


if ~isempty(indx)
    for kk = 1:length(indx)
        kernel = ones(1,nkernel(indx(kk)));
        kernel = kernel./norm(kernel,1);
        idx = cM(kk):cM(kk+1)-1;
        data = abs(x(idx)).^2;
        W(idx) = sqrt(conv(data,kernel,'same'));
    end
end
        
if strcmp(param.method,'pew')
    sol = pew(x,lambda,W);
else
    sol = st(x,lambda,W);
end
