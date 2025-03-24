function plot_predicted_vs_observed(obj, varargin)
% plot_predicted_vs_observed Plot predicted versus observed outcomes.
%
% For regression (more than two unique values of Y):
%   - Plots a scatter plot of predicted (.yfit) vs. observed (.Y) values.
%   - Adds a regression line.
%   - Calculates the Pearson correlation and p-value.
%   - Prints a statement whether the regression is significant (p < 0.05).
%
% For classification (exactly two unique values of Y):
%   - Creates a cell array with predicted values separated by true class.
%   - Calls barplot_columns (from CANlab core tools) to generate violin plots.
%
% :Usage:
% ::
%     pm_obj.plot_predicted_vs_observed();
%     pm_obj.plot_predicted_vs_observed('noplot');
%     pm_obj.plot_predicted_vs_observed('color', [0.2 0.4 0.8]);
%     pm_obj.plot_predicted_vs_observed('color', {[0.2 0.4 0.8] [.8 .4 .2]});
%
% :Optional Inputs:
%
%   **'noplot':**
%        If provided, suppresses all graphical output.
%
%   **'color':**
%        An [r g b] triplet specifying the plot color for regression scatter.
%        Default is medium-blue.
%        OR, for classification models, a cell array with two [r g b]
%        triplets, one for each class
%
% :Outputs:
%
% :Examples:
% ::
%     % Regression example:
%     pm_obj.predicted_observed_scatterplot();
%     % Classification example:
%     pm_obj.predicted_observed_scatterplot('noplot');
%
% (Additional plotting and analysis methods are forthcoming.)

% Default settings:
doPlot = true;
plotColor = [0 0.4470 0.7410];  % medium-blue, for regression

% for classification, use default colors from seaborn_colors (assumed available from CANlab core tools).
colors = seaborn_colors(10);
colors = colors([4 8]);

% Parse optional arguments.
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'noplot'
                doPlot = false;
            case {'color', 'colors'}
                if i < length(varargin)
                    plotColor = varargin{i+1};
                end
        end
    end
end

% Check that true outcomes and predictions exist.
if isempty(obj.Y) || isempty(obj.yfit)
    error('predictive_model:MissingData', 'Both .Y and .yfit must be defined.');
end

uniqueY = unique(obj.Y);

if numel(uniqueY) > 2
    %% Regression Case
    % Compute correlation and p-value.
    [r, p] = corr(obj.Y, obj.yfit, 'Rows','complete');
    if p < 0.05
        sigStr = 'significant';
    else
        sigStr = 'not significant';
    end
    fprintf('Regression: r = %.3f, p = %.3f. The regression is %s.\n', r, p, sigStr);

    if doPlot
        N = numel(obj.Y);
        % Set point size inversely proportional to sqrt(N) (with a minimum size).
        pointSize = max(20, 100/sqrt(N));

        scatter(obj.yfit, obj.Y, pointSize, plotColor, 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        hold on;
        % Fit and plot a linear regression line.
        pFit = polyfit(obj.Y, obj.yfit, 1);
        xFit = linspace(min(obj.Y), max(obj.Y), 100);
        yFit = polyval(pFit, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 2);
        hold off;

        % Label axes using .Y_name and .X_name if they exist.
        if isprop(obj, 'Y_name') && ~isempty(obj.Y_name)
            ylabel(obj.Y_name);
        else
            ylabel('Observed Outcome');
        end
        if isprop(obj, 'X_name') && ~isempty(obj.X_name)
            xlabel(obj.X_name);
        else
            xlabel('Predicted Outcome (Model scores)');
        end
        title('Predicted vs. Observed Scatterplot');
    end

else
    %% Classification Case (Two unique outcomes)
    % Create a cell array with predictions grouped by true outcome.
    cellData = cell(2,1);
    cellData{1} = obj.dist_from_hyperplane_xval(obj.Y == uniqueY(1));
    cellData{2} = obj.dist_from_hyperplane_xval(obj.Y == uniqueY(2));

    if doPlot

        if length(obj.dist_from_hyperplane_xval) > 500
            extra_args = {'noind'};
        end

        % barplot_columns is assumed to plot violin and bar plots.
        cla
        barplot_columns(cellData, 'colors', colors, 'nobars', 'nofig', 'names', obj.class_labels, extra_args{:});

        xlabel('True class')

        if isprop(obj, 'X_name') && ~isempty(obj.X_name)
            ylabel(obj.X_name);
        else
            ylabel('Model scores');
        end

        title('Predicted Outcomes by True Class');

        if isprop(obj, 'class_labels') && ~isempty(obj.class_labels)

            if length(obj.class_labels) ~= 2
                error('class_labels must be cell array with two cells, each containing a string')
            end

            set(gca, 'XTickLabel', obj.class_labels, 'XTickLabelRotation', 45);

        end

        if isprop(obj, 'Y_name') && ~isempty(obj.Y_name)

            xlabel(obj.Y_name);
        else
            xlabel('True class');
        end
        title('Predicted Outcomes by True Class');

        % Print string with t-test of class differences

        if obj.mult_obs_within_person
            disp('t-test for diffs in model scores, unpaired observations (single-interval)')
            disp('Not defined yet. Add for paired observations and update me!!')

        else
            disp('t-test for diffs in model scores, unpaired observations (single-interval)')
            ttest2_printout(cellData{1}, cellData{2});

        end

    end % doplot


    % fprintf('Classification problem detected (two unique outcomes). Barplot of predictions generated.\n');
end

set(gcf, 'Color', 'w')
set(gca, 'FontSize', 18)

end % function
