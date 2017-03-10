function [x,L] = pad(x,L)
% :Usage:
% ::
%
%     x = pad(x,L)
% 
% pads x with zeros of length L
%
% or, if L is a vector, to longest of two (NOT DONE)
%
% COLUMN VECTORS

if length(L) == 1
    
    x = [x; zeros(L,size(x,2))];
    
else
    len = max(size(x,1),size(L,1));
    x = [x; zeros(len - size(x,1),size(x,2))];
    L = [L; zeros(len - size(L,1),size(L,2))];
    
end

return
