function cm = confusionchart(obj, varargin)
% confusionchart  Confusion chart of cross-validated predictions.
%
% Thin wrapper around MATLAB's confusionchart using the model's true
% outcomes (obj.Y) and cross-validated predictions
% (obj.fitted_values.yfit). Class tick labels come from obj.class_labels
% when set.
%
% :Usage:
% ::
%     cm = confusionchart(pm);
%     cm = confusionchart(pm, 'Normalization', 'row-normalized');
%
% :Inputs:
%
%   **obj:**
%        a cross-validated @predictive_model (classification task).
%
% :Optional Inputs (name/value, forwarded to MATLAB confusionchart):
%
%   Any options accepted by MATLAB's confusionchart, e.g. 'Normalization'
%   ('absolute' | 'row-normalized' | 'column-normalized'), 'Title'.
%
% :Outputs:
%
%   **cm:**
%        the ConfusionMatrixChart handle.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, dat.dat', dat.Y);
%     confusionchart(pm, 'Normalization', 'row-normalized');
%
% :See also:
%   confusion_matrix, plot, rocplot, predictive_model

    if ~isstruct(obj.fitted_values) || ~isfield(obj.fitted_values, 'yfit') ...
            || isempty(obj.fitted_values.yfit) || isempty(obj.Y)
        error('predictive_model:confusionchart:MissingData', ...
            'Need obj.Y and obj.fitted_values.yfit. Run crossval() (or fit()) first.');
    end

    Y    = obj.Y(:);
    yfit = obj.fitted_values.yfit(:);

    % Default to row-normalized unless the caller passed 'Normalization'.
    args = varargin;
    if ~any(cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'Normalization'), args))
        args = [args, {'Normalization', 'row-normalized'}];
    end

    cm = confusionchart(Y, yfit, args{:});

    if ~isempty(obj.class_labels) && numel(obj.class_labels) == numel(unique(Y(~isnan(Y))))
        try, cm.ClassLabels = obj.class_labels; end %#ok<NOCOM,TRYNC>
    end
    cm.OffDiagonalColor = [1 1 1];
    cm.DiagonalColor    = [.2 .5 1];
    cm.FontSize         = 16;
    set(gcf, 'Color', 'w');
end
