function [c,g,shift,M,Ls,fb,tgtfl ] = mywvltanalysis( sig, fmin, Fs, bins, bw, win )
%function needed for audioinpainting to produce vector ouput.

[ccc,g,shift,M,Ls,fb,tgtfl] = wvlttrans(sig,fmin,Fs,bins,bw,win);
c = c2mat(ccc,M);


end

