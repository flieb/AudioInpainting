function [g,shift,M,fb] = nswvltwp2win(fmin,bw,bins,fs,L)

%function to build frame for the wp2inp 

%desired bandwidth in samples
bwd = bw*L/fs;
%desired fmin position in samples
posbg = fmin*L/fs;

%get first window with desired bandwidth bwd
[g1,BL,BR] = getwin(bwd,bins);
[~,maxg] = max(g1);

%shift parameter for first window such that the maximum is located at the
%desired fmin position
pos = round(posbg - maxg);

%compute the corresponding left and right bounds
BR = BR + pos;
BL = BL + pos;

%compute number of scales, BR should not be greater than L/2
scales = floor(log2(L/2/BR)*bins)+1;

%compute powers of 2
pow2 = 2.^((0:scales-1)/bins)';

%compute bandwidths of all scaled windows as well as left bound of all
%windows
bwds = bwd*pow2;
blv = floor(BL.*pow2);


%empty container for results
g = cell(2*scales + 2,1);
M = zeros(2*scales+2,1);
shift = zeros(2*scales + 2,1);


%number of shifts is equal to bandwidths
M(2:scales+1) = ceil(bwds);
M(end:-1:scales+3) = M(2:scales+1);


% get windows with desired bandwidths bwds
[g(2:scales+1),tmp_bl] = arrayfun( @(x) getwin(x,bins),bwds,'UniformOutput',0);
%scale them:
g(2:scales+1) = arrayfun( @(x) g{x}/sqrt(M(x)), 2:scales+1,'UniformOutput',0);
%reverse them for upper half:
g(end:-1:scales+3) = cellfun(@(x) flipud(x),g(2:scales+1),'UniformOutput',0);
%get length of all window functions (note that the length is not equal to
%the bandwidth
Lg = cellfun(@length,g(2:scales+1));
%do ifftshift
g = cellfun(@(x) ifftshift(x),g,'UniformOutput',0);


%compute shift parameter
shift(2:scales+1) = blv - cell2mat(tmp_bl) + floor(Lg./2) + 2;
shift(end:-1:scales + 3) = L - shift(2:scales+1) + 2 + mod(Lg-1,2);
shift(1) = 1;
shift(scales + 2) = floor(L/2) + 1;
shift = [L-mod(shift(end)-1,L); diff(shift)];
posit = cumsum(shift)-shift(1);


% sum diagonal of frame operator
diagonal = zeros(L,1);
for ii = 2:scales+1
  win_range = mod(posit(ii)+...
      (-floor(length(g{ii})/2):ceil(length(g{ii})/2)-1)-1,L)+1;
  diagonal(win_range) = diagonal(win_range) + M(ii)*fftshift(g{ii}).^2 ;   
end
for ii = scales+3:2*scales+2
  win_range = mod(posit(ii)+...
      (-floor(length(g{ii})/2):ceil(length(g{ii})/2)-1)-1,L)+1;
  diagonal(win_range) = diagonal(win_range) + M(ii)*fftshift(g{ii}).^2 ;
end


%set number of points to cover 0 and Nyquist frequencies
LowLim = 2*BR;
UppLim = L/2 - blv(end);

%get corresponding min and max 
maxXL = max(diagonal(1:LowLim));
maxXU = max(diagonal(blv(end):floor(L/2)));


% Get window for 0 Hz
M(1) = 2*LowLim;
g{1} = sqrt(max((maxXL-diagonal([1:LowLim,end-LowLim+1:end]))/M(1),0));

%get window for Nyquist frequency
M(scales+2) = 2*UppLim;
g{scales+2} = sqrt(max((maxXU-diagonal(floor(L/2)+[1:UppLim,-UppLim+1:0]))/M(scales+2),0));


% Compute the frame bounds if output parameter 'fb' is desired
if nargout == 4

    % Finalize the construction of the diagonal 
    Lg = cellfun(@length,g);
    for ii = [1,scales+2]
      range = mod(posit(ii)+(-floor(Lg(ii)/2):ceil(Lg(ii)/2)-1)-1,L)+1;
      diagonal(range) = diagonal(range) + (fftshift(g{ii}).^2)*M(ii);
    end
    fb = [min(diagonal) max(diagonal)];
end





 
 
 



