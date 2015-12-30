function [y,ntrimmed,allw] = trimts(y,sd,X,varargin)
% 1.  Adjusts for scan effects (unless X is empty)
%
% 2.  Windsorizes timeseries to sd standard deviations
%       - Recursive: 3 steps
%
% 3.  Adds scan effects back into timeseries
%
% :Usage:
% ::
%
%     function [y,ntrimmed,allw] = trimts(y,sd,X,[do spike correct],[trimming iterations],[MADs])
%
% Spike correct: Some attempt at automatic adjustment for abrupt level
% shifts in data; default is 0, enter 1 to do this
% "Spikes" are IDd as values more than 10 MADs (by default) from moving average with 20 image FWHM
% Replaces "spike" data values with moving average
%
% iterations: number of cycles through trimming; default is 3
% 
% MADs: allows you to change the number of MADs above which values are IDd
% as "spikes"
%
% ..
%    Modified 2/9/05, 02/2008
%    Tor Wager
% ..
% filter y using X matrix; yf is residuals

ntrimmed = 0;
allw = [];
MADs = 10;

if nargin < 3, X = []; end

if ~isempty(X),
    mfit = X * pinv(X) * y;
    yf = y - mfit;
else
    yf = y;
end

niter = 3; spike = 0;
if ~isempty(varargin), spike = varargin{1}; end
if length(varargin) > 1, niter = varargin{2}; end
if length(varargin) > 2, MADs = varargin{3}; end

if spike
    
    % Figure out which values deviate from a moving average
    basefit = moving_average('gaussian', yf, 20);
    
    devs = abs(yf - basefit);
    
    % More than 10 MADs from moving average with 20 image FWHM
    allw = find( devs ./ median(devs) > MADs ) ;
    
    yf(allw)  = basefit(allw);
    
%     create_figure('ts'); plot(yf); hold on; plot(basefit,'r');
%     plot_vertical_line(find(wh));
    
% % %     
% % %     % attempt to correct for session-to-session baseline diffs
% % %     
% % %     tmp = diff(yf);
% % %     mad12 = median(abs(tmp)) * Inf;    % robust est of change with time (dt)
% % %     wh = find(abs(tmp) > mad12);
% % %     n = 20;
% % %     
% % %     for i = 1:length(wh), 
% % %         st = max(wh(i) - (n-1),1);  % start value for avg
% % %         en = max(wh(i),1);
% % %         st2 = wh(i)+1;
% % %         en2 = min(wh(i)+n,length(yf));  % end value for after
% % %         wh2 = st2:en2;
% % %         m = mean(yf(wh(i)+1:en2)) - mean(yf(st:en)); % average of 5 tp after -  5 time points before
% % %         %figure;plot(st:en2,yf(st:en2));
% % %         yf(wh(i)+1:end) = yf(wh(i)+1:end) - m;,
% % %     end
% % %     
% % %     
% % %     % do spike correction!  Interpolate values linearly with 1 nearest
% % %     % neighbor
% % %     %
% % %     % replace first val with mean
% % %     n = min(length(yf),50);
% % %     yf(1) = mean(yf(1:n));
% % %     
% % %     tmp = diff(yf);
% % %     mad5 = median(abs(tmp)) * 5;    % robust est of tail of dist. of change with time (dt)
% % %     wh = find(abs(tmp) > mad5);
% % %     
% % %     % find paired changes that are w/i 3 scans
% % %     whd = diff(wh);
% % %     wh = wh(whd < 3);
% % %     whd = whd(whd < 3);
% % % 
% % %     % value of spike is avg of pre-spike and post-spike val.
% % %     wh(wh == 1 | wh == length(yf)) = [];
% % %     for i = 1:length(wh)
% % %         if wh(i)+1+whd(i)<=length(yf);
% % %         yf(wh(i)+1) = mean([yf(wh(i)) yf(wh(i)+1+whd(i))]);
% % %         %else
% % %        % yf(wh(i)+1) = mean([yf(wh(i)) yf(length(yf))]);    
% % %     end
% % %     end
    

end
    
    
    
% trim residuals to sd standard deviations
% "Windsorize"

my = mean(yf);

for i = 1:niter
    
    if all(yf - nanmean(yf) < eps), return, end     % no variance in timeseries!!
    
    yf2 = scale(yf);
    w = find(abs(yf2) > sd);
    yf(w) = my + sd * std(yf) * sign(yf(w)-my);
    
    allw = [allw; w];
end

allw = unique(allw);

% put means back into yf
if ~isempty(X)
    y = yf + mfit;
else
    y = yf;
end

ntrimmed = length(allw);  % w) + length(w2);

return
