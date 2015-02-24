function Wcomp = contrast2d(W,comp)
%function Wcomp = contrast2d(W,comp)
%
% W is a 2D matrix of states nested within participants
% comp is a column vector with contrast weights to apply to each
% participant
%
% e.g.,
% contrast2d(W,[1 -1]')
% contrast2d(DATA.INDSCAL.W,DATA.SPEC.comps(1,:)')


comp = repmat(comp,1,size(W,2));

tmp = [];
for i = 1:size(comp,1):size(W,1)
    tmp(end+1,:) = sum(W(i:i+size(comp,1)-1,:) .* comp);
end

Wcomp = tmp;

return