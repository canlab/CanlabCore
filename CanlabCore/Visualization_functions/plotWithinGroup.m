function plotWithinGroup(dataGroups, groupLabels, varargin)
% Plot multiple groups of data with options for connecting lines, half-violin plots, scatter plots, and box plots.
% Michael Sun, Ph.D.
% 
% USAGE:
%   plotWithinGroup(dataGroups, groupLabels, varargin)
% 
% INPUTS:
%   dataGroups: Cell array, each cell contains data for one group
%   groupLabels: Cell array of strings, labels for each group
%   varargin: Optional parameters to control plot appearance and components
% 
% OUTPUT:
%   None
% 
% EXAMPLES:
%   % Example 1: Basic usage
%   dataGroup1 = randn(100, 1);
%   dataGroup2 = randn(100, 1) + 1;
%   dataGroups = {dataGroup1, dataGroup2};
%   groupLabels = {'Group 1', 'Group 2'};
%   plotWithinGroup(dataGroups, groupLabels);
% 
%   % Example 2: Customized appearance
%   plotWithinGroup(dataGroups, groupLabels, 'PlotScatter', false, 'PlotBox', false);
% 
%   % Example 3: Specify colors
%   customColors = [0.1 0.2 0.3; 0.4 0.5 0.6]; % RGB values
%   plotWithinGroup(dataGroups, groupLabels, 'Colors', customColors);
% 
%   % Example 4: Adjust offset
%   plotWithinGroup(dataGroups, groupLabels, 'Offset', 0.2);
% 
%   % Example 5: Change density scale factor
%   plotWithinGroup(dataGroups, groupLabels, 'DensityScaleFactor', 0.01);
% 
%   % Example 6: Disable specific components
%   plotWithinGroup(dataGroups, groupLabels, 'PlotLines', false, 'PlotHalfViolin', false);
% 
%   % Example 7: Change mid-violin style
%   plotWithinGroup(dataGroups, groupLabels, 'MidViolin', 'face2face');
%
% NOTES:
%   This function plots multiple groups of data with various options for visualization. It can be useful for comparing distributions of data across different groups.

% Validate inputs
nGroups = length(dataGroups);
assert(nGroups == length(groupLabels), 'Number of data groups must match number of labels');

% Default settings
defaultColors = lines(nGroups); % MATLAB's 'lines' colormap for up to 7 groups, then repeats
offset = 0.15;
densityScaleFactor = 0.005;
plotLines = true;
plotHalfViolin = true;
plotScatter = true;
plotBox = true;
midViolin = 'fullviolin'; % The other option is 'face2face'

% Parse optional inputs
p = inputParser;
addParameter(p, 'Colors', defaultColors, @(x) ismatrix(x) && size(x,2) == 3);
addParameter(p, 'IndColors', defaultColors, @(x) ismatrix(x) && size(x,2) == 3);
addParameter(p, 'Offset', offset, @isscalar);
addParameter(p, 'DensityScaleFactor', densityScaleFactor, @isscalar);
addParameter(p, 'PlotLines', plotLines, @islogical);
addParameter(p, 'PlotHalfViolin', plotHalfViolin, @islogical);
addParameter(p, 'PlotScatter', plotScatter, @islogical);
addParameter(p, 'PlotBox', plotBox, @islogical);
addParameter(p, 'MidViolin', midViolin, @ischar);
parse(p, varargin{:});

% Extract parsed results
colors = p.Results.Colors;
indcolors = p.Results.IndColors;
offset = p.Results.Offset;
densityScaleFactor = p.Results.DensityScaleFactor;
plotLines = p.Results.PlotLines;
plotHalfViolin = p.Results.PlotHalfViolin;
plotScatter = p.Results.PlotScatter;
plotBox = p.Results.PlotBox;
midViolin = p.Results.MidViolin;

% Create figure
% figure; 
% 
hold on;

% Plot each group
for i = 1:nGroups
    currentData = dataGroups{i};
    timePoint = i;
    
    % Plot connecting lines between groups if enabled
    if plotLines && i > 1
        for j = 1:length(dataGroups{i})
            prevData = dataGroups{i-1};
            % line([i-1, i] + [offset, -offset], [prevData(j), currentData(j)], 'Color', colors(i,:), 'LineWidth', 1);
            line([i-1, i] + [offset, -offset], [prevData(j), currentData(j)], 'Color', indcolors(j,:), 'LineWidth', 1);
        end
    end
    
    % Plot half-violin if enabled
    if plotHalfViolin
        [pdf, value] = ksdensity(currentData);
        pdf = (pdf / max(pdf)) * max(abs(currentData) - min(currentData)) * densityScaleFactor;
        if i == 1
            xCoords = [timePoint * ones(1, length(value)), timePoint - pdf(end:-1:1)];
            fill(xCoords, [value, value(end:-1:1)], colors(i,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
        elseif i == nGroups % Adjust for last group if needed
            xCoords = [timePoint * ones(1, length(value)), timePoint + pdf(end:-1:1)];
            fill(xCoords, [value, value(end:-1:1)], colors(i,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
        else
            if strcmpi(midViolin, 'fullViolin')
                xCoords = [timePoint * ones(1, length(value)), timePoint - pdf(end:-1:1)];
                fill(xCoords, [value, value(end:-1:1)], colors(i,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
                xCoords = [timePoint * ones(1, length(value)), timePoint + pdf(end:-1:1)];
                fill(xCoords, [value, value(end:-1:1)], colors(i,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
            end
            if strcmpi(midViolin, 'face2face')
                xCoords = [timePoint * ones(1, length(value)), timePoint - pdf(end:-1:1)];
                fill(xCoords+(offset/1.5), [value, value(end:-1:1)], colors(i,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
                xCoords = [timePoint * ones(1, length(value)), timePoint + pdf(end:-1:1)];
                fill(xCoords-(offset/1.5), [value, value(end:-1:1)], colors(i,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
            end

        end
        
    end
    
    % Plot scatter if enabled
    if plotScatter

    
        if i ~= 1 && i ~= nGroups
            % Plot middle-groups
            scatter(timePoint + offset, currentData, 36, indcolors, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'jitter','on', 'jitterAmount',0.04);
            scatter(timePoint - offset, currentData, 36, indcolors, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'jitter','on', 'jitterAmount',0.04);
        else
            % Plot edge-groups
            scatter(timePoint + (i==1)*offset - (i==nGroups)*offset, currentData, 36, indcolors, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1, 'jitter','on', 'jitterAmount',0.04);
        end


    end
    
    % Plot box plot if enabled
    if plotBox
        boxplot(currentData, 'Positions', timePoint, 'Widths', 0.1, 'Colors', colors(i,:), 'Symbol', '', 'Whisker', 1);
        % Ensure the boxplots are in the background
        uistack(findobj(gca,'Type','boxplot'), 'bottom');
    end


end

% Final plot adjustments
xlim([0.5, nGroups + 0.5]);
xticks(1:nGroups);
xticklabels(groupLabels);
ylabel('Data Values');
title('Multi-group Data Plot');
hold off;
end
