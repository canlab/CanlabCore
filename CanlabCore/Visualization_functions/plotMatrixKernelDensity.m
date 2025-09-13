function axis_handles = plotMatrixKernelDensity(vals, minPlotVal, maxPlotVal, colors, names)
% axis_handles = plotMatrixKernelDensity plots the KDE for each column of vals in subplots.
%
%   plotMatrixKernelDensity(vals, minPlotVal, maxPlotVal) takes a matrix where
%   each column is a dataset and calls kernelDensitySmoothedHistogram for each column,
%   placing the plots in a single figure with one subplot (row) per column.
%
%   Optional:
%       colors - a cell array of 1x3 RGB vectors (one per column) specifying the line
%                color for each plot. If provided, each subplot will use the corresponding
%                color for both the KDE line and the shaded area.
%
% Example:
%   plotMatrixKernelDensity(dataMatrix, 0, 10, {[0.9 0.3 0.6], [0.1 0.7 0.2], [0.4 0.4 0.8]}, {'ref_map1' 'ref_map2' 'ref_map3'});

    % Validate input.
    if ~ismatrix(vals)
        error('Input "vals" must be a matrix.');
    end
    
    [~, numCols] = size(vals);
    
    % If colors are not provided, use default
    if nargin < 4
        colors = seaborn_colors(size(vals, 2));
    end

    % If names are not provided, use default
    if nargin < 5
        names = {}; for i = 1:numCols, names{i} = ''; end
    end

    names = format_strings_for_legend(names);

    % Create a new figure.
    create_figure('matrix_density_histograms');
    
    % Loop over each column and plot in its own subplot.
    for i = 1:numCols
        % Create a subplot: one row per column.
        ax = subplot(numCols, 1, i);
        axis_handles(i) = ax;
        
        % If a color is provided for this column, pass it along.
        if ~isempty(colors) && numel(colors) >= i && ~isempty(colors{i})

            kernelDensitySmoothedHistogram(vals(:, i), minPlotVal, maxPlotVal, ...
                'color', colors{i}, 'shade', true, 'shadecolor', colors{i}, ...
                'axishandle', ax);
        else
            kernelDensitySmoothedHistogram(vals(:, i), minPlotVal, maxPlotVal, ...
                'axishandle', ax);
        end
        
        % Remove extra white space around the axes.
        set(ax, 'LooseInset', get(ax, 'TightInset'));

        if i < numCols 
            xlabel(''); 
            set(gca, 'XColor', 'none')
        end

        title(names{i})
        set(gca, 'XLim', [minPlotVal maxPlotVal])

        % formatting for larger matrices
        if numCols > 8
            title('')
            ylabel('')
            set(gca, 'YTickLabel', names{i}, 'YTick', 0)
        end

    end

end % main function



