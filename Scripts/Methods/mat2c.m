function [c] = mat2c(vec_c,M)

% inverse to c2mat

    ind = 1;
    c = cell(length(M),1);

    for ii = 1:length(M)
        c{ii} = vec_c(ind : ind + M(ii)-1);
        ind = ind + M(ii);
    end


end