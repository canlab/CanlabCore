function [y2,ntrimmed,spikes, yy] = splinetrim(y,varargin)
% function [y,ntrimmed,spikes, yfit] = splinetrim(y,[iqrmult],[knotrate],[X],['p'])
%
% Uses a robust measure of deviations in a timeseries gradient
% to find high-velocity 'spikes', presumed to be artifacts
%
% Uses spline interpolation to replace spikes with reasonable values.
%
% y         a timeseries
%
% optional inputs:
% iqrmult   how many times the interquartile range above which velocities are
%           outliers, default is 1.5
% knotrate  sets knot points every k observations, default is 3
% X         matrix of session means or other linear regressors to remove 
% 'p'       plot the results
% X and 'p' can be entered in any order, but after iqrmult and knotrate
%
% Example (good for eye tracking):
% [y2,nt] = splinetrim(trialdat,3,5,'p'); nt
%
%  3/29/05, Tor Wager



% -------------------------------------------------
% * setup
% -------------------------------------------------

iqrmult = 3;  % how many times the interquartile range above which pts are outliers
knotrate = 5;   % in images
X = [];
plotme = 0;


if length(varargin) > 0, iqrmult = varargin{1};, end
if length(varargin) > 1, knotrate = varargin{2};, end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
        case 'p', plotme = 1;
        otherwise error('Unrecognized argument: acceptable is ''p'' for plot.')
        end
    else
        if ndims(varargin{i}) == 2 & sum(size(varargin{i}))>2
            X = varargin{i};
        end
    end
end
% filter y using X matrix; yf is residuals
if ~isempty(X),
    mfit = X * pinv(X) * y;
    y = y - mfit;
end


% ---------------------------------------------------------
% find spike regions
% ---------------------------------------------------------
% EXCLUDE FROM CONSIDERING AS KNOT POINTS BASED ON VELOCITY
veloc = gradient(y); veloc(end) = 0;    % clamp last value to 0
tmp = abs(veloc - median(veloc));

% threshold is mad + 1.5 times the interquartile range of abs. deviations
thr = median(tmp) + iqrmult * (prctile(tmp,75) - prctile(tmp,25));

spikes = tmp>thr;


% EXCLUDE FROM CONSIDERING AS KNOT POINTS BASED ON OUTLIER STATUS
d = abs(y - median(y));    % distance from mean
thr = median(d) + iqrmult * (prctile(d,75) - prctile(d,25));
sptmp = d>thr;
spikes = spikes + sptmp;


% smooth this some, so that low-deviation regions in the middle of spikes
% do not get counted
spikes2 = smooth_timeseries(spikes,4);  % excluded from being knot points

spikes = find(spikes>0);        % which points to interpolate in the end
spikes2 = find(spikes2>0);      % larger set of pts to not include as knot pts


% -------------------------------------------------
% * get spline fit 
% -------------------------------------------------
bp = zeros(size(y));
bp(1:knotrate:length(y)) = 1;  
bp(isnan(y)) = 0;                % ignore NaNs  
%bp(end) = 1;                     % clamp last point to be an endpt
bp(spikes2) = 0;                 % do not put knot points on spikes

bp = find(bp);

nbp = length(bp);
bpy = zeros(nbp, 1);

% figure out medians for each segment - 
% this is the interp knot point.
ytmp = y;
ytmp(spikes2) = NaN;  % get rid of spikes

for i = 1:nbp
    
    if i == 1
        st = 1;
    else
        st = bp(i) - round((bp(i) - bp(i-1)) ./ 2);
    end
    
    if i == nbp
        en = length(y);
    else
        en = bp(i) + round((bp(i+1) - bp(i)) ./ 2);
    end
    data = ytmp(st:en);
    bpy(i) = median(data);

end


% xx = 1:length(y);
yy = spline(bp, bpy, 1:length(y));
if ~iscol(yy), yy = yy'; end


% -------------------------------------------------
% * find pts that are really far from spline fit
% -------------------------------------------------
d = abs(y - yy);    % distance from spline fit
thr = median(d) + iqrmult * (prctile(d,75) - prctile(d,25));

spikes3 = d>thr;

spikes = unique([spikes; find(spikes3>0)]);        % final answer
% replace y data with spline fits where spikes occur
y2 = y;
y2(spikes) = yy(spikes);



ntrimmed = length(spikes); 


if plotme
    figure; hold on;
    plot(y,'k','LineWidth',2);
    plot(y2,'g');
    plot(yy,'b');
    plot(bp, yy(bp), 'b.','MarkerSize', 10);
    legend({'Original' 'Adjusted' 'Spline fit'})
    plot(spikes,y2(spikes),'ro','MarkerFaceColor','r','MarkerSize',6);
end

return