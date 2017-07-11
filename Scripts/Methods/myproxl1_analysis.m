function sol = myproxl1_analysis(x, lambda, param)

% proximal operator for l1 norm => softthresholding
%analysis version of the inverse problem:

%softthresholding:
st = @(s, T) max(abs(s)-T,0)./(max(abs(s)-T,0)+T).*s;
ew = @(s,T) s.*max(1-(T./abs(s)).^2,0);
%st = @(s,T) s.*max(1-T./abs(s),0);
%hard thresholding:
ht = @(s, T) (abs(s)>T).*s;

if nargin < 3
    param.verbose = 1;
end
if ~isfield(param,'verbose'), param.verbose = 1; end


tmp = param.A(x);
sol = x + param.At(st(tmp,lambda) - tmp);




%sol = st(x,lambda);
%sol = ew(x,lambda);


if param.verbose
    fprintf(' prox_L1: ||x||_1 = %e\n', norm(x,1));
end


