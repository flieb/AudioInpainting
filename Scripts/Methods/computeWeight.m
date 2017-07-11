function W = computeWeight(x, param)

%% compute the weights

%this only works, as long as the groups are singular in one direction
[xr,xc] = size(x);
[r,c] = size(param.group); %r -> frequency groups, c -> time groups
if param.overlap == 0
    if r==1 %groups in time direction
        %M = mean(reshape(abs(x)',c,numel(x)/c));
        M = sqrt(sum( reshape(abs(x)',c,numel(x)/c).^2,1 ));
        W = reshape(repmat(M,c,1),xc,xr)';
    elseif c==1 %groups in freq direction
        M = mean(reshape(abs(x),r,numel(x)/r));
        W = reshape(repmat(M,r,1),xr,xc);
    end
else

    %W = convFreqWeights(x.',r,c).';
    neigh = ones(r,c);
    neigh = neigh./(norm(neigh(:),1));
    windowcenter = ceil(r/2);
    windowcenter2= ceil(c/2);

    W2 = zeros(xr + r-1,xc+c-1,'double');
    %extend the borders
    W2(windowcenter:xr+windowcenter-1,windowcenter2:xc+windowcenter2-1)=abs(x).^2;
    W2(1:windowcenter-1,:) = flipud( W2(windowcenter:2*(windowcenter-1),:) );
    W2(xr+windowcenter:end,:) = flipud( W2(xr-r+2*windowcenter :xr+windowcenter-1,:) );
    W2(:,1:windowcenter2-1)= fliplr( W2(:,windowcenter2:2*(windowcenter2-1)));
    W2(:,xc+windowcenter2:end)= fliplr( W2(:,xc-c+2*windowcenter2:xc+windowcenter2-1));

    %do convolution
    W = (conv2(W2,neigh,'valid')).^(1/2);
end