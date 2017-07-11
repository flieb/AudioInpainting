function W = compWeights(X,param)

%computes the Weights for group lasso
%X are the tf coefficients
%G is a matrix of same size with integers, each integer belongs to one
%group

%% compute the weights

G = param.groupmat;
size_g = size(G);
G = G(:);
W = zeros(size(G));
numgroups = max(G(:));

if ~isfield(param, 'w')
    param.w = ones(1,numgroups);
end


for kk=1:param.numgroups
   
    W(param.groupidx{kk}) = sqrt(param.w(kk)).*sqrt(sum(X(param.groupidx{kk}).^2));
end
   
W = reshape(W,size_g);
