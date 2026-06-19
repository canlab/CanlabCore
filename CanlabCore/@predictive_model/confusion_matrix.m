function [rawConf, normConf] = confusion_matrix(obj, varargin)
% confusion_matrix Compute and (optionally) plot the confusion matrix.
%
% :Usage:
% ::
%     [rawConf, normConf] = obj.confusion_matrix();
%     [rawConf, normConf] = obj.confusion_matrix('noplot');
%
% :Inputs:
%
%   **'noplot':** (Optional) String flag that, if provided, suppresses the plot.
%
% :Outputs:
%
%   **rawConf:** Confusion matrix of raw counts.
%   **normConf:** Confusion matrix normalized to percentages of the number of
%                observations in each true class.
%
% This method computes the confusion matrix using the object's true outcomes (Y)
% and predicted outcomes (yfit). It uses MATLAB's built-in confusionmat function
% to compute the raw counts and then normalizes each row to percentages.
%
% If no 'noplot' flag is provided, the method also displays a confusion chart using
% MATLAB's confusionchart function.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, dat.dat', dat.Y);
%     [raw, norm] = confusion_matrix(pm);            % compute and plot
%     [raw, norm] = confusion_matrix(pm, 'noplot');  % compute only
%
% :See also:
%   confusionchart, plot, rocplot

    % Default: plot the confusion matrix.
    plotFlag = true;

    for i = 1:length(varargin)
        if ischar(varargin{i}) && strcmpi(varargin{i}, 'noplot')
            plotFlag = false;
        end
    end

    % Fetch predictions from the categorised sub-struct (the flat .yfit
    % alias was removed during property deduplication).
    if isstruct(obj.fitted_values) && isfield(obj.fitted_values, 'yfit')
        yfit = obj.fitted_values.yfit;
    else
        yfit = [];
    end

    % Ensure true outcomes (Y) and predictions (yfit) are defined.
    if isempty(obj.Y) || isempty(yfit)
        error('predictive_model:MissingData', ...
            'Both true outcomes (Y) and predicted outcomes (fitted_values.yfit) must be defined.');
    end

    % Compute the confusion matrix using MATLAB's built-in function.
    rawConf = confusionmat(obj.Y, yfit);

    % Normalize each row to percentages (per true class).
    normConf = rawConf;
    for i = 1:size(rawConf, 1)
        rowSum = sum(rawConf(i, :));
        if rowSum > 0
            normConf(i, :) = (rawConf(i, :) / rowSum) * 100;
        else
            normConf(i, :) = 0;
        end
    end

    % If plotting is not suppressed, display the confusion chart.
    if plotFlag
        create_figure('confchart');
        hold off;

        cm = confusionchart(obj.Y, yfit, 'Normalization', 'row-normalized');  % doesn't work for me... , 'classLabels', {'TrueClass','FalseClass'}); % , 'ClassLabels', ClassLabels);
        cm.OffDiagonalColor = [1 1 1];
        cm.DiagonalColor = [.2 .5 1];
        cm.FontSize = 18;

    end
end