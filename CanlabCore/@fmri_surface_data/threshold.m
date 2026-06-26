function obj = threshold(obj, thresh, varargin)
% threshold Raw-value threshold of grayordinate data (zeros sub-threshold values).
%
% :Usage:
% ::
%     obj = threshold(obj, [lo hi])
%     obj = threshold(obj, t, 'positive')
%     obj = threshold(obj, t, 'negative')
%
% Surface analogue of the raw-value branch of image thresholding: grayordinate
% values inside the threshold range are set to 0; .dat keeps its full size (D5b).
%
% Cluster-extent thresholding is available via the 'k' option (mesh-graph
% connected components for cortex, 26-connectivity for subcortex; see
% reparse_contiguous).
%
% :Inputs:
%   **obj:**    an fmri_surface_data object.
%   **thresh:** scalar t (keeps |value| >= t, two-tailed) or [lo hi] (keeps
%               values <= lo OR >= hi; i.e. removes the open interval (lo,hi)).
%
% :Optional Inputs:
%   **'positive':** keep only values >= t (with scalar t).
%   **'negative':** keep only values <= -abs(t) (with scalar t).
%   **'k', N:**     cluster-extent threshold: after the raw-value threshold,
%                   remove contiguous clusters smaller than N grayordinates
%                   (applied per map column).
%
% :Outputs:
%   **obj:** thresholded object (sub-threshold grayordinates zeroed).
%
% :See also: fmri_surface_data, apply_mask, reparse_contiguous

direction = 'two';
if any(strcmpi(varargin, 'positive')), direction = 'pos'; end
if any(strcmpi(varargin, 'negative')), direction = 'neg'; end
kextent = 0;
for i = 1:numel(varargin)
    if (ischar(varargin{i}) || isstring(varargin{i})) && strcmpi(varargin{i}, 'k') && i < numel(varargin)
        kextent = varargin{i+1};
    end
end

d = double(obj.dat);

if numel(thresh) == 2
    lo = thresh(1); hi = thresh(2);
    keep = d <= lo | d >= hi;
else
    t = abs(thresh);
    switch direction
        case 'pos', keep = d >= t;
        case 'neg', keep = d <= -t;
        otherwise,  keep = abs(d) >= t;
    end
end

d(~keep) = 0;
obj.dat = single(d);

% Cluster-extent: drop clusters smaller than k grayordinates, per map
if kextent > 1
    for col = 1:size(obj.dat, 2)
        tmp = reparse_contiguous(obj, 'which_image', col);
        cl = tmp.brain_model.cluster;
        if isempty(cl), continue; end
        szs = accumarray(cl(cl > 0), 1, [max(cl) 1]);
        small = find(szs < kextent);
        drop = ismember(cl, small);
        dd = obj.dat(:, col); dd(drop) = 0; obj.dat(:, col) = dd;
    end
    obj.history{end+1} = sprintf('threshold cluster-extent k>=%d applied', kextent);
end

obj.history{end+1} = sprintf('threshold (raw-value): kept %.1f%% of values', 100*mean(keep(:)));
end
