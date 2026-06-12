function varargout = plot(obj, varargin)
% plot  Default visualization for a predictive_model (task-dispatched).
%
% Draws the core performance figure for the model's task and, when the model
% carries them, additional figures for the artifacts you have computed:
%
%   Core (always):
%     regression      -> predicted-vs-observed scatter with regression line
%                        (plot_predicted_vs_observed)
%     classification  -> two panels: model scores by true class (violin)
%                        and the cross-validated ROC curve
%
%   Conditional (only if present, unless suppressed):
%     weight map      -> montage(pm) in a new figure, when a weight map has
%                        been cached (pm.weights.weight_obj, e.g. via
%                        weight_map_object / fmri_data.predict 'newapi').
%                        Add 'surface' to also render the cortical surface.
%     permutation     -> plot_permutation(pm) in a new figure, when
%                        permutation_test has populated pm.permutation_results.
%
% This is a convenience entry point. The underlying methods
% (plot_predicted_vs_observed, rocplot, confusionchart, montage, surface,
% plot_permutation) can be called directly for finer control.
%
% :Usage:
% ::
%     plot(pm);
%     plot(pm, 'noplot');                 % compute/report only
%     plot(pm, 'color', [.2 .4 .8]);      % forwarded to the scatter/violin
%     plot(pm, 'noweights');              % skip the weight-map montage
%     plot(pm, 'noperm');                 % skip the permutation histogram
%     plot(pm, 'surface');                % also render the weight surface
%
% :Inputs:
%
%   **obj:**
%        a fitted / cross-validated @predictive_model.
%
% :Optional Inputs:
%
%   **'noweights':**
%        do not draw the weight-map montage even if a weight map is cached.
%
%   **'surface':**
%        also draw the cortical-surface rendering of the weight map.
%
%   **'noperm':**
%        do not draw the permutation null histogram even if available.
%
%   Other options ('noplot', 'color', ...) are forwarded to
%   plot_predicted_vs_observed.
%
% :Outputs:
%
%   **h:**
%        (optional) handle of the core performance figure.
%
% :Examples:
% ::
%     % Classification, with weights + permutation auto-shown
%     dat = load_image_set('DPSP_hotwarm', 'noverbose');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, X, Y, 'groups', id);
%     pm  = weight_map_object(pm, dat);            % cache a weight map
%     pm  = permutation_test(pm, X, Y, 'nperm', 1000, 'groups', id);
%     plot(pm);                                    % ROC + violin, montage, perm hist
%
%     % Regression: predicted-vs-observed scatter
%     pm = predictive_model('algorithm','svr','task','regression');
%     pm = crossval(pm, X, Y);
%     plot(pm);
%
% :See also:
%   plot_predicted_vs_observed, rocplot, confusionchart, montage, surface,
%   plot_permutation

    % --- pull out plot-control flags this method handles itself ---
    [varargin, do_weights] = pop_flag(varargin, 'noweights', true);   % default show
    [varargin, do_surface] = pop_flag(varargin, 'surface',   false);  % default off
    [varargin, do_perm]    = pop_flag(varargin, 'noperm',    true);   % default show

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

    % ---- core performance figure ----
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

    if ~do_plot
        if nargout > 0, varargout{1} = h; end
        return
    end

    % ---- conditional: weight-map montage (+ optional surface) ----
    has_weight_obj = isstruct(obj.weights) && isfield(obj.weights, 'weight_obj') ...
        && ~isempty(obj.weights.weight_obj);
    if do_weights && has_weight_obj
        try
            create_figure('predictive_model weights'); axis off
            montage(obj);
        catch ME
            fprintf('plot: weight montage skipped (%s)\n', ME.message);
        end
        if do_surface
            try
                create_figure('predictive_model weight surface'); axis off
                surface(obj);
            catch ME
                fprintf('plot: weight surface skipped (%s)\n', ME.message);
            end
        end
    end

    % ---- conditional: permutation null histogram ----
    has_perm = isstruct(obj.permutation_results) ...
        && isfield(obj.permutation_results, 'null_scores') ...
        && ~isempty(obj.permutation_results.null_scores);
    if do_perm && has_perm
        try
            create_figure('predictive_model permutation');
            plot_permutation(obj);
        catch ME
            fprintf('plot: permutation histogram skipped (%s)\n', ME.message);
        end
    end

    if nargout > 0, varargout{1} = h; end
end


% -------------------------------------------------------------------------
function [args, val] = pop_flag(args, name, default)
% Remove a bare flag from args; val = ~default if present, else default.
    hit = cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, name), args);
    if any(hit)
        args(hit) = [];
        val = ~default;
    else
        val = default;
    end
end
