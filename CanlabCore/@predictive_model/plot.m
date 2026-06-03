function varargout = plot(obj, varargin)
% plot  Default visualization for a predictive_model (task-dispatched).
%
% Dispatches on obj.task:
%   regression      -> predicted-vs-observed scatter with regression line
%                      (plot_predicted_vs_observed)
%   classification  -> a two-panel figure: model scores by true class
%                      (violin) and the cross-validated ROC curve.
%
% This is a convenience entry point. The underlying methods
% (plot_predicted_vs_observed, rocplot, confusionchart, montage, surface)
% can be called directly for finer control.
%
% :Usage:
% ::
%     plot(pm);
%     plot(pm, 'noplot');                 % compute/report only
%     plot(pm, 'color', [.2 .4 .8]);      % forwarded to the scatter/violin
%
% :Inputs:
%
%   **obj:**
%        a fitted / cross-validated @predictive_model.
%
% :Optional Inputs:
%
%   Forwarded to plot_predicted_vs_observed (e.g. 'noplot', 'color').
%
% :Outputs:
%
%   **h:**
%        (optional) the figure handle.
%
% :Examples:
% ::
%     % Classification
%     dat = load_image_set('DPSP_hotwarm');
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, dat.dat', dat.Y);
%     plot(pm);
%
%     % Regression
%     pm = predictive_model('algorithm','svr','task','regression');
%     pm = crossval(pm, X, Y);
%     plot(pm);
%
% :See also:
%   plot_predicted_vs_observed, rocplot, confusionchart, montage, surface

    % Determine task (fall back to Y cardinality).
    task = obj.task;
    if isempty(task)
        if ~isempty(obj.Y) && numel(unique(obj.Y(~isnan(obj.Y)))) <= 2
            task = 'classification';
        else
            task = 'regression';
        end
    end

    do_plot = ~any(cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'noplot'), varargin));

    if strcmpi(task, 'classification') && do_plot
        h = create_figure('predictive_model plot', 1, 2);

        subplot(1, 2, 1);
        plot_predicted_vs_observed(obj, varargin{:});

        subplot(1, 2, 2);
        try
            rocplot(obj);   % roc_plot draws into the current axes (gca)
        catch ME
            title('ROC unavailable');
            fprintf('plot: ROC panel skipped (%s)\n', ME.message);
        end
    else
        % Regression, or 'noplot' (let plot_predicted_vs_observed report).
        h = gcf;
        plot_predicted_vs_observed(obj, varargin{:});
    end

    if nargout > 0, varargout{1} = h; end
end
