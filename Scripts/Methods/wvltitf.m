function f = wvltitf(c,gd,am)

%inverse wavelet transform using dual window
%L = size(gd,2);

if ( sum(am) == size(gd,1) )
    %full inverse transform
    f = real(ifft(sum(bsxfun(@times,fft(c,[],2),gd))));
else
    if ( size(gd,1) > size(gd,2) ) %means gd is transposed for the mex file
        f = real(ifft(comp_wvltitf2(cellfun(@fft, c, 'UniformOutput',false),gd,am)))';
    else %regular case
        f = zeros(1,size(gd,2));
        for kk=1:size(gd,1)
            f = f + gd(kk,:).*repmat(fft(c{kk}),am(kk),1).';
        end
        f = real(ifft(f))';
    end

end

