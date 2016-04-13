function [line_handle patch_handle] = plot_error(varargin)
% Plot a matrix with shaded error area around it
%
% :Usage:
% ::
%
%     [line_handle patch_handle] = plot_error(varargin)
%
%
% :Examples:
% :: 
%
%     % plots matrix Y against x-values in X, where Y is a matrix with 
%     % each row representing a signal. The shaded area represents the 
%     % standard error across the columns of Y.
%     PLOT_ERROR(X, Y) 
%
%     % plots the data matrix Y versus its index.
%     PLOT_ERROR(Y) 
%
%     % uses external data for the error areas. In this case, Y is assumed 
%     % to be a mean timeseries already. Y and errorData must be vectors 
%     % of the same length.
%     PLOT_ERROR(..., 'errorData', errorData)
%
%     % The following example indicates whether to handle NaNs in the data. 
%     % If set, plot_error will use nanmean, nanstd, etc. Off by default.
%     PLOT_ERROR(..., 'allowNaNs', [0|1]) 
%
%     % The following example plots the line according to the designated 
%     % ColorSpec string, and shades the error area by the color of the line
%     PLOT_ERROR(..., colorSpecString) 
%
%     % The following example plots into the axes designated by the AX 
%     % axes handle.
%     PLOT_ERROR(AX, ...) 
%
%     % The following example returns the handle of the main line object
%     [line_handle patch_handle] = PLOT_ERROR 

    using_external_error_data = 0;
    plotArgs = {};
    meanFun = @mean;
    stdFun = @std;
    skip_next_arg_for_plot = 0;

    if(isempty(varargin))
        error('No arguments to %s', mfilename);
    end

    args = varargin;
    if(ishandle(args{1}))
        h = args{1};
        if(~strcmp(get(h, 'Type'), 'axes'))
            error('Handle %d must be an axes handle.', h);
        end
        args(1) = [];
    end

    if(isempty(args{1}) || (~isvector(args{1}) && ~ismatrix(args{1})))
        error('First argument (if not an axes handle) must be a matrix or vector');
    elseif(length(args) > 1 && isvector(args{1}) && ismatrix(args{2}) && length(args{1}) == size(args{2}, 2))
        x = args{1};
        y = args{2};
        args(1:2) = [];
    else
        y = args{1};
        x = 1:size(y, 2);
        args(1) = [];
    end
    x = x(:);

    for i=1:length(args)
        if(ischar(args{i}))
            switch(args{i})
                case 'errorData'
                    if(~isvector(y) || ~isvector(args{i+1}) || (length(y) ~= length(args{i+1})))
                        error('If passing in external error data, Y and the error data must be vectors of equal length');
                    end
                    meanData = y;
                    errorData = args{i+1};
                    using_external_error_data = 1;
                    skip_next_arg_for_plot = 1;
                case 'allowNaNs'
                    meanFun = @nanmean;
                    stdFun = @nanstd;
                    skip_next_arg_for_plot = 1;
                otherwise
                    if(skip_next_arg_for_plot)
                        skip_next_arg_for_plot = 0;
                    else
                        plotArgs{end+1} = args{i};
                    end
            end
        else
            if(skip_next_arg_for_plot)
                skip_next_arg_for_plot = 0;
            else
                plotArgs{end+1} = args{i};
            end
        end
    end

    if(~using_external_error_data)
        if(isvector(y))
            error('Y is a vector. Nothing meaningful can be drawn. Either pass in a matrix, or pass in external error data to plot around Y.');
        end
        meanData = meanFun(y);
        errorData = stdFun(y) / sqrt(size(y, 1));
    end

    if(~exist('h', 'var') || isempty(h))
        h = gca();
    end
    line_handle = plot(h, x, meanData, plotArgs{:});
    hold(h, 'on');
    xVertices = [x; flipud(x)];
    yVertices = [meanData+errorData meanData(end:-1:1)-errorData(end:-1:1)];
    patch_handle = fill(xVertices, yVertices, get(line_handle, 'Color'), 'EdgeColor', 'none', 'FaceAlpha', .25, 'Parent', h);
    hold(h, 'off');
end


% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%
% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 && c>1)
    ismat = 1;
  else
    ismat = 0;
  end

end
