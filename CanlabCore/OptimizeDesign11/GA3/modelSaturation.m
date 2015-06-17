function model = modelSaturation(model,varargin)
% function model = modelSaturation(model,[threshold],[soft])
%
% replaces all values above some threshold with a specified value
% to model saturation effects
% default = 4 (4 x max HRF)
%
% a 3rd argument signals a soft thresholding, where
% the argument is the exponent (ex) in a thresholding model:
% y(y > thresh) = thresh - (thresh .^ ex) + (y(y > thresh) .^ ex);
%
% 3/30/01 Tor Wager
%

if nargin > 1, thresh = varargin{1};,else thresh = 4;,end

if nargin < 3
    model(model(:,:) > thresh) = thresh;
else
    if isempty(varargin{2}), ex = .5;, else, ex = varargin{2};,end
    model(model>thresh) = 2 - (2 .^ ex) + (model(model>thresh) .^ ex);
end

return