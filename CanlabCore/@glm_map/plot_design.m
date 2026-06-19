function plot_design(obj, varargin)
% Plot the design of a glm_map object.
%
% If the model is an event/1st-level design with a small number of event
% types of interest (<= 12), the events of interest are shown as color-coded
% regressor line plots with boxes marking each event and its duration (via
% plotDesign), alongside a heat map of the full design matrix (including any
% nuisance covariates and the intercept). Otherwise (direct/group designs, or
% many event types) only the full-design heat map is shown. Heat maps use
% 'YDir','reverse' (observation 0 at the top) and axis tight.
%
% :Usage:
% ::
%
%     plot_design(obj)
%     plot_design(obj, 'max_event_types', 12)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design available (obj.X non-empty).
%
% :Optional Inputs:
%
%   **'max_event_types', [k]:**
%        Use the plotDesign line-plot format only when the number of event
%        types of interest is <= k (default 12).
%
% :See also:
%   - plotDesign, diagnostics, fmri_glm_design_matrix.plot
%
% ..
%    2026 - Hooks into plotDesign for events of interest; full-design heat map.
% ..

X = obj.X;
if isempty(X)
    error('glm_map:NoDesign', 'No design matrix available (obj.X is empty).');
end

max_event_types = 12;
wh = find(strcmpi(varargin, 'max_event_types'));
if ~isempty(wh), max_event_types = varargin{wh(1) + 1}; end

% -------------------------------------------------------------------------
% Can we draw the event-of-interest line plot? (event design, onsets present,
% and a manageable number of event types of interest)
% -------------------------------------------------------------------------
[ons, names] = local_interest_onsets(obj);
do_plotdesign = ~isempty(ons) && numel(ons) <= max_event_types;

if do_plotdesign
    create_figure('glm_map design', 1, 2);

    subplot(1, 2, 1);
    % Durations are carried as the 2nd column of each onset cell (seconds).
    plotDesign(ons, [], obj.TR, 'samefig');
    title(sprintf('Events of interest (%d types)', numel(ons)));
    xlabel('Time (s)');
    if ~isempty(names) && numel(names) <= 12
        legend(names, 'Interpreter', 'none', 'Location', 'best');
    end

    subplot(1, 2, 2);
    local_heatmap(obj, X);
else
    create_figure('glm_map design');
    local_heatmap(obj, X);
end

end % plot_design


% =========================================================================
% Local helpers
% =========================================================================
function [ons, names] = local_interest_onsets(obj)
% Collect onsets (sec, with durations as a 2nd column) and names for the
% conditions of interest of session 1. Returns {} if not an event design.
[ons, names] = deal({});

if isempty(obj.design) || ~isa(obj.design, 'fmri_glm_design_matrix'), return, end
if isempty(obj.design.Sess) || ~isfield(obj.design.Sess(1), 'U') || isempty(obj.design.Sess(1).U), return, end

U = obj.design.Sess(1).U;
TR = obj.TR;

% onset/duration units -> seconds
to_sec = 1;
if ~isempty(obj.design.xBF) && isfield(obj.design.xBF(1), 'UNITS') ...
        && any(strcmpi(obj.design.xBF(1).UNITS, {'scans', 'tr', 'trs'}))
    to_sec = TR;
end

for i = 1:numel(U)
    o = U(i).ons(:) * to_sec;

    dd = [];
    if isfield(U(i), 'dur'), dd = U(i).dur; end
    if isempty(dd), dd = 0; end
    dd = dd(:);
    if isscalar(dd), dd = repmat(dd, numel(o), 1); end
    dd = dd * to_sec;

    ons{i}   = [o dd];                     % onsets + durations (sec)
    nm = U(i).name; if iscell(nm), if isempty(nm), nm = sprintf('Cond%d', i); else, nm = nm{1}; end, end
    names{i} = nm;
end
end


function local_heatmap(obj, X)
% Heat map of the full design matrix (incl. nuisance covariates / intercept).
imagesc(X);
colormap(gray);
colorbar;
set(gca, 'YDir', 'reverse');               % observation 0 at the top
axis tight;
xlabel('Regressor');
ylabel('Image / observation');
title(sprintf('Full design matrix (%d interest, %d nuisance, %d intercept)', ...
    sum(obj.wh_interest), sum(obj.wh_nuisance), sum(obj.wh_intercept)));

rn = obj.regressor_names;
if ~isempty(rn) && numel(rn) == size(X, 2)
    % Mark nuisance (n) and intercept (i) columns in the tick labels
    labels = rn(:)';
    for j = 1:numel(labels)
        if obj.wh_intercept(j),    labels{j} = [labels{j} ' (icpt)'];
        elseif obj.wh_nuisance(j), labels{j} = [labels{j} ' (nuis)'];
        end
    end
    set(gca, 'XTick', 1:size(X, 2), 'XTickLabel', labels, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
end
end
