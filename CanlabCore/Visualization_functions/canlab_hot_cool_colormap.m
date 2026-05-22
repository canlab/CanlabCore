function cmap = canlab_hot_cool_colormap(varargin)
% Hot/cool diverging colormap with white at zero.
%
% Returns a colormap that goes from deep blue/purple at the most negative
% end, fades to pure white at zero, and rises to orange/yellow at the most
% positive end. Designed to be paired with a symmetric color axis
% (CLim = [-m m]) so the value 0 always maps to white. All-positive data
% then naturally renders in white-to-yellow, all-negative data in
% blue-to-white, and mixed-sign data uses the full diverging range.
%
% :Usage:
% ::
%
%     cmap = canlab_hot_cool_colormap
%     cmap = canlab_hot_cool_colormap(n)
%
% :Optional Inputs:
%
%   **n:**
%        Number of colormap levels. Must be even so that exactly half the
%        levels are below and half above zero. Default: 256.
%
% :Outputs:
%
%   **cmap:**
%        n-by-3 RGB colormap. Rows 1..n/2 cover negative values (deep blue
%        → light blue → white), rows n/2+1..n cover positive values
%        (white → orange → yellow).
%
% :Examples:
% ::
%
%     % Apply to current axes (pair with symmetric CLim)
%     m = max(abs(get(gca, 'CLim')));
%     set(gca, 'CLim', [-m m]);
%     colormap(canlab_hot_cool_colormap);
%
%     % Use as the default for a display_slices figure
%     load_image_set('emotionreg');
%     display_slices(mean(ans));
%     colormap(canlab_hot_cool_colormap);
%
% :See also:
%   - display_slices
%   - colormap
%
% ..
%    Tor Wager, May 2026
% ..

if isempty(varargin) || isempty(varargin{1})
    n = 256;
else
    n = varargin{1};
end

if mod(n, 2) ~= 0
    n = n + 1;  % force even so half is below 0 and half above
end

half = n / 2;

% Negative half: deep blue/purple -> pale -> white
neg_anchors = [0.00 0.00 0.50;   % deep blue (most negative)
               0.10 0.15 0.85;   % blue
               0.35 0.20 0.85;   % blue-purple
               0.70 0.65 0.95;   % pale lavender
               1.00 1.00 1.00];  % white at 0

% Positive half: white -> orange -> yellow
pos_anchors = [1.00 1.00 1.00;   % white at 0
               1.00 0.90 0.55;   % pale cream
               1.00 0.65 0.10;   % orange
               1.00 0.40 0.00;   % deep orange
               1.00 1.00 0.00];  % yellow (most positive)

neg = interp1(linspace(0, 1, size(neg_anchors, 1)), neg_anchors, ...
              linspace(0, 1, half), 'linear');
pos = interp1(linspace(0, 1, size(pos_anchors, 1)), pos_anchors, ...
              linspace(0, 1, half), 'linear');

cmap = [neg; pos];
cmap = max(0, min(1, cmap));

end
