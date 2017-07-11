function [f] = mywvltsynthesis( c_vec,g,shift,M,Ls,fb,tgtfl )
%function needed for audioinpainting

c = mat2c(c_vec,M);
[f] = invwvlttrans(c,g,shift,M,tgtfl,fb,Ls);



end