function sol = myproxl21(x, lambda, param)


% proximal operator for l21 norm => group lasso with non overlapping groups

if nargin < 3
    param.verbose = 1;
    error('no group structure specified, using groups of size 1!');
end
if ~isfield(param,'verbose'), param.verbose = 1; end

st = @(s,T,W) s.*max(1-T./W,0);
pew = @(s,T,W) s.*max(1-(T./W).^2,0);

%% compute the weights

W = computeWeight(x,param);

%%
if strcmp(param.method,'pew')
    sol = pew(x,lambda,W);
else
    sol = st(x,lambda,W);
end
%sol = pew(x,lambda,W);
% figure(3), subplot(131); plotdgtreal(x,100,2000,44100,60,'nocolorbar');
% subplot(132); plotdgtreal(W,100,2000,44100,60,'nocolorbar');
% subplot(133); plotdgtreal(sol,100,2000,44100,60,'nocolorbar');
% pause();

if param.verbose
    fprintf(' prox_L1: ||x||_21 = %e\n', norm(W,1));
end