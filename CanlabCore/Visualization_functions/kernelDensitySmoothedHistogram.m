function [density_vals, xValues] = kernelDensitySmoothedHistogram(vals, minPlotVal, maxPlotVal, varargin)
% kernelDensitySmoothedHistogram plots a kernel density-smoothed histogram.
%
%   kernelDensitySmoothedHistogram(vals, minPlotVal, maxPlotVal) plots the
%   kernel density estimate (KDE) of the data in vector vals, draws a horizontal
%   gray line from minPlotVal to maxPlotVal at y = 0, and uses a default line color.
%
%   Optional name-value pair arguments:
%       'bandwidth'  - A scalar controlling the amount of smoothing (bandwidth)
%                      used in the KDE. Default: [] (automatic selection).
%       'color'  - A 1x3 vector specifying the RGB color for the KDE line.
%                      Default: [0.9 0.3 0.6].
%       'shade'      - A logical flag indicating whether to shade the area under
%                      the KDE curve. Default: true.
%       'shadecolor' - A 1x3 vector specifying the RGB color for the shaded area.
%                      Default: same as color.
%       'shadealpha' - A scalar between 0 and 1 for the transparency of the shaded
%                      area. Default: 0.5.
%       'axishandle'  - An axes handle in which to plot the histogram. If not provided,
%                       a new figure is created.
%
% Example:
%   kernelDensitySmoothedHistogram(data, 0, 10, 'bandwidth', 0.5, ...
%       'color', [0.1 0.7 0.2], 'shade', true);

    % Validate that vals is a vector.
    if ~isvector(vals)
        error('Input "vals" must be a vector.');
    end

    if nargin < 2
        
        minPlotVal = min(vals);
        maxPlotVal = max(vals);

    end

    % Parse optional parameters.
    p = inputParser;
    addParameter(p, 'bandwidth', [], @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'color', [0.9, 0.3, 0.6], @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'shade', true, @(x) islogical(x) || (isnumeric(x) && (x==0 || x==1)));
    addParameter(p, 'shadecolor', [], @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'shadealpha', 0.5, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    addParameter(p, 'axishandle', [], @(x) isempty(x) || ishandle(x));
    parse(p, varargin{:});
    
    bandwidth = p.Results.bandwidth;
    color = p.Results.color;
    shadeFlag = logical(p.Results.shade);
    shadealpha = p.Results.shadealpha;
    axishandle = p.Results.axishandle;

    % If shadecolor is not provided, default to color.
    if isempty(p.Results.shadecolor)
        shadecolor = color;
    else
        shadecolor = p.Results.shadecolor;
    end

    % Determine the number of evaluation points (bins) using a heuristic.
    nBins = max(10, round(sqrt(numel(vals))));
    
    % Compute the kernel density estimate (KDE).
    if ~isempty(bandwidth)
        [density_vals, xValues] = ksdensity(vals, 'NumPoints', nBins, 'bandwidth', bandwidth);
    else
        [density_vals, xValues] = ksdensity(vals, 'NumPoints', nBins);
    end

    % Ensure column vectors for consistency.
    xValues = xValues(:);
    density_vals = density_vals(:);

    % Create a new figure and hold for multiple plots, if we're not using
    % existing axis handle
    if isempty(axishandle)
        figure;
        hold on;
    else
        axes(axishandle)
        hold on;
    end


    % Optionally add shading under the KDE curve.
    if shadeFlag
        % Create patch coordinates: the curve and the baseline (y=0).
        patch([xValues; flipud(xValues)], [density_vals; zeros(size(density_vals))], shadecolor, ...
              'FaceAlpha', shadealpha, 'EdgeColor', 'none');
    end

    % Plot the KDE as a line.
    plot(xValues, density_vals, 'Color', color, 'LineWidth', 3);

    % Draw a horizontal gray line from minPlotVal to maxPlotVal at y = 0.
    plot([minPlotVal, maxPlotVal], [0, 0], 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2);

    hold off;
    xlabel('Value');
    ylabel('Density');
    title('Kernel Density-Smoothed Histogram');
    set(gca, 'FontSize', 18)

end % Main function



