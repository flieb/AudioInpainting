function Wd = nsgl(data,kernelsize)

%nonstationary group lasso with no overlap. 
%helper function for nsgl. handles one scale

l = length(data);
m = mod(l,kernelsize);

data2 = data(1:end-m);


%version1
W = sum(reshape(data2',kernelsize,floor(l/kernelsize)),1);
Wd = reshape(repmat(W,kernelsize,1),size(data2,1),size(data2,2));

if m
    tmp = sum(data(end-m+1:end));
    Wd(end+1:end+m) = repmat(tmp,1,m);
end

