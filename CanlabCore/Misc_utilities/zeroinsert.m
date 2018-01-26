function yout = zeroinsert(wasbad, y)
% Re-insert removed CASES (rows) and fill with zeros
%
% :Usage:
% ::
%
%     yout = zeroinsert(wasbad, y)
%
% wasbad is indicator for removed cases, of size size(yout).  y is data.
%
% if you enter y', inserts VARIABLES (cols).  here, pass in v x n matrix, y'
% to fill empty/removed vars
%
% See nanremove.m and naninsert.m

wh = find(wasbad);

if isempty(wh), yout = y; return, end

if length(wh) + size(y, 1) ~= length(wasbad)
    disp('Illegal removed-cases vector. length must equal size(y, 1) + # left-out cases.');
    fprintf('Left out = %3.0f, size(y, 1) = %3.0f, sum = %3.0f, length(wasbad) = %3.0f', length(wh), size(y, 1), length(wh)+size(y, 1), length(wasbad));
    error('Quitting')
end
    
yout = zeros(length(wasbad), size(y, 2));

yout(~wasbad, :) = y;

% yout(1:wh(1) - 1, :) = y(1:wh(1) - 1, :);
% 
% for i = 2:length(wh)
%     ystart = wh(i-1) + 2 - i; % NaN index value - num previous removed; wh(i-1) + 1 - i + 1;
%     yend = wh(i) - i; % wh(i) - 1 - i + 1
%     
%     if yend > size(y, 1)
%         error('Illegal value of %3.0f data index, which exceeds data rows of %3.0f. Bad removed-cases vector?', yend, size(y, 1));
%     end
%         
%     yout(wh(i-1) + 1 : wh(i) - 1, :) = y(ystart : yend, :);
% end
% % last segment
% i = i+1;
% if isempty(i), i = 2; end
% ystart = wh(end) + 2 - i;
% 
% if ystart <= size(y, 1)
%     yout(wh(end)+1:end, :) = y(ystart:size(y, 1), :);
% end

end
