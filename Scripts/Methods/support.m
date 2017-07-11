function [ il, ir, supp] = support( win, tol )
% approximates support of function win where 1 - tol percent of the norm
% are located. 
% i_l: left index
% i_r: right index
% support: i_r - i_l + 1

[r,c] = size(win);
if r>c
    win = win.';
end

n = length(win);


if nargin < 2
    tol = 1e-6;
end

if norm(win) ~= 1
    win = win / norm(win);
end


win2 = abs(win).^2;

tmp = cumsum(win2);
il = find(tmp >= tol/2,1,'first');
tmp2 = cumsum(fliplr(win2));
ir = n-(find(tmp2 >= tol/2,1,'first')-1);

supp = ir-il+1;

% i_left = 1;
% i_right = 1;
% normleft = sum(win2);
% normright = normleft;
% temp = normleft - tol;
% 
% while(normleft > temp)
%     normleft = sum(win2(i_left : end));
%     i_left = i_left + 1;
% end
% 
% while(normright > temp)
%     normright = sum(win2(1 : end - i_right));
%     i_right = i_right + 1;
% end
% 
% i_l = i_left - 2; %why -2 : -1 then normleft < 1- 1e-6, and -1 such that normleft > 1- 1e-6
% i_r = n - (i_right - 2);
%
%
% supp = i_r - i_l + 1;


end

