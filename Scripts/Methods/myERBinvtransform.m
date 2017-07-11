function [f] = myERBinvtransform( c_vec,gdERB,shiftERB,MERB,ls)
%function needed for audioinpainting



cERB = mat2c(c_vec,MERB);

[f] = nsigtf(cERB,gdERB,shiftERB,ls);



end