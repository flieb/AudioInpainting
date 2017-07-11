function [vec_c] = c2mat(c,M)

% function to rearange cell output from NSGT toolbox
% to a one-dimensional vector.

    n = sum(M);
    ind = 1;
    T = zeros(n,1);

    for ii = 1:length(c)
        T(ind : ind + M(ii) - 1) = c{ii};
        ind = ind + M(ii);
    end

    vec_c = T;
end