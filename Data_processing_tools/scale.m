function xout = scale(x,varargin)
% x = scale(x,[just center])
%
% centers and scales column vectors
% to mean 0 and st. deviation 1

%Programmer's notes:
%   2/17/14 - Luke: rewrote using vectorization instead of loop -MUCH FASTER!
%                   Still removes nans like original function
%   6/13/14 - Wani: after Luke's edit, mean-centering didn't work right.
%                   Fixed it. 

xout = NaN .* zeros(size(x));

[nanvec x_no_nan] = nanremove(x);  %removes an entire row if any nan

% no data.  return original input.
if isempty(x_no_nan)
    xout = x;
    return
end

x_no_nan = x_no_nan - repmat(mean(x_no_nan),size(x_no_nan,1),1); %remove mean

if isempty(varargin) || varargin{1} == 0
    x_no_nan = bsxfun(@rdivide,x_no_nan,max(eps,std(x_no_nan))); %divide by std
end

xout(~nanvec,:) = x_no_nan;


%TOR'S OLD SCALE USING LOOP - replaced by LC: 2/27/14
% xout = NaN .* zeros(size(x));
%
% for i = 1:size(x,2)
%
%     [nanvec x_no_nan] = nanremove(x(:,i));
%
%     if isempty(x_no_nan)
%         % no data.  return original input.
%         xout = x;
%         return
%     end
%
%     x_no_nan = x_no_nan - mean(x_no_nan); %repmat(mean(x_no_nan),size(x_no_nan,1),1);
%
%     if length(varargin) == 0 || varargin{1} == 0
%
%         x_no_nan = x_no_nan ./ max(eps, std(x_no_nan)); %repmat(std(x_no_nan),size(x_no_nan,1),1);
%
%     end
%
%     xout(~nanvec,i) = x_no_nan;
% end

return