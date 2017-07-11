function c = wvlttf(f,g,am)

%wavelet transform according to my algorithm...

[c,r] = size(f);
if c>r, f = f.'; end


%full transform: means am not needed...
%tic
%c = ifft( bsxfun(@times,g,fft(f)) ,[],2);
%toc

%g = g(1:end/2+1,:);
if ( sum(am) == size(g,1) )
    %full transform
    c = ifft( bsxfun(@times,g,fft(f)) ,[],2);
else
    fh = fft(f);
    M = size(g,2)./am;
    
    if ( size(g,1) > size(g,2) )
        c = cellfun(@ifft,comp_wvlttf(fh,g,am),'UniformOutput',0);
    else
        c = arrayfun(@(x) ifft(sum(reshape(g(x,:).*fh,M(x),am(x)),2))./am(x),1:size(g,1),'UniformOutput',0);
    end
    
    
    
%     tic
%     c = cell(size(g,1),1); 
%     for kk = 1:size(g,1)
%         c{kk} = (sum(reshape(g(kk,:).*fh,M(kk),am(kk)),2))./am(kk);
%     end
%     c = cellfun(@ifft,c,'UniformOutput',0);
%     toc
    
    
    %compnorm(c{100},c2{100});
end