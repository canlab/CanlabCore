function [h, z] = plot_joint_hist_contour(z, xbins, ybins, color, varargin)
% plot a 95% 2-D density region for a 2-D histogram
%
% :Usage:
% ::
%
%    h = plot_joint_hist_contour(z, xbins, ybins, color, ['confval', confval], ['maxalpha', maxalpha])
%
% :Inputs:
%
%   **z:**
%        2-D histogram values, counts in bins.  See joint_hist.m
%
%   **confval:**
%        optional input; [0 - 1], retain this proportion of values in confidence region
%
%   **maxalpha:**
%        optional; [0 - 1], maximum transparency 
%
% :Outputs:
%
%   **h:**
%        handle to graphical contour object
%
%   **z:**
%        thresholded z matrix of counts
%
% :Examples:
% ::
%
%    z = joint_hist(nnmfscores{i}{j}(:, 1),nnmfscores{i}{j}(:, 2), 50, 'noplot');
%    h = plot_joint_hist_contour(z, [0 0 1]);

confval = .95;
maxalpha = 1;

%if ~isempty(varargin), confval = varargin{1}; end

for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                    % functional commands
                case {'confval', 'alpha'}, confval = varargin{i+1};
                case 'maxalpha', maxalpha = varargin{i+1};
                
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

clear cnt
for i = 1:max(z(:))
    cnt(i) = sum(z(z <= i));
end

crit = find(cnt <= (1 - confval) * nansum(z(:)));

%%
z(z <= max(crit)) = NaN;

[~, h] = contourf(xbins, ybins, z);

set(h, 'LineColor', 'none')

hh = get(h, 'Children');

%%
vals = get(hh, 'CData');
vals = cat(1, vals{:});
vals = maxalpha .* vals ./ max(vals);

%set(hh, 'FaceColor', 'b')
%set(hh, 'FaceAlpha', 'interp')

for i = 1:length(hh)
    
    if ~isnan(vals(i))
        set(hh(i), 'FaceColor', color, 'FaceAlpha', vals(i))
    end
    %set(hh(i), 'AlphaDataMapping', 'scaled', 'AlphaData', get(hh(i), 'CData'), 'FaceAlpha', 'interp')
    %set(hh(i), 'AlphaDataMapping', 'scaled', 'AlphaData', get(hh(i), 'CData'))
    %set(hh(i), 'FaceAlpha', 'interp')
end

end
