function yout = naninsert(nanvec, y)
% Insert NaNs back into data series from which they were removed
%
% :Usage:
% ::
%
%     yout = naninsert(nanvec, y)
%
% nanvec is logical indicator for NaNs, of size size(yout).  y is data.
%
% If y is a matrix, length(nanvec) should equal number of rows in y.
% NaNs are inserted for true values in nanvec in all columns of y.
%
% ..
%    Tor Wager, July 2007; Bug fix, March 2008;
%    Updated functionality to include all columns in matrix, April 2015
%  ..
%
% :See also:
% nanremove.m

    wh = find(nanvec);
    
    if isempty(wh), yout = y; return, end
    
    yout = NaN * ones(length(nanvec), size(y, 2));

    yout(~nanvec, :) = y;
    
%     yout(1:wh(1) - 1) = y(1:wh(1) - 1);
%     for i = 2:length(wh)
%         ystart = wh(i-1) + 2 - i; % NaN index value - num previous removed; wh(i-1) + 1 - i + 1; 
%         yend = wh(i) - i; % wh(i) - 1 - i + 1
%         
%         yout(wh(i-1) + 1 : wh(i) - 1, :) = y(ystart : yend, :);
%     end
%     % last segment
%         i = i+1;
%         if isempty(i), i = 2; end
%         ystart = wh(end) + 2 - i;
%         yout(wh(end)+1:end, :) = y(ystart:size(y, 1), :);
%     
%     yout(wh, :) = NaN;

end  % function
