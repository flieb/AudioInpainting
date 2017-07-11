function [c] = myERBlettransform( sig,gERB,shiftERB,MERB)
%function needed for audioinpainting to produce vector ouput.
%Input:
%[gERB,shiftERB,MERB] = nsgerbwin(bins,Fs,ls,'Qvar',Qvar);
%gdERB = nsdual(gERB,shiftERB,MERB);

cERB = nsgtf(sig,gERB,shiftERB,MERB);
c = c2mat(cERB,MERB);


end
