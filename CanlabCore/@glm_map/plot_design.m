function plot_design(obj, varargin)
% Plot the design of a glm_map object.
%
% If the model is an event/1st-level design with a small number of event
% types of interest (<= 12), each event type is drawn in its own panel showing
% the object's *actual* basis-convolved regressor(s) -- one line per basis
% function, so the real basis set is respected -- with boxes marking each event
% and its duration. A heat map of the full design matrix (including any
% nuisance covariates and the intercept) is shown alongside. Otherwise
% (direct/group designs, or many event types) only the full-design heat map is
% shown. Heat maps use 'YDir','reverse' (observation 0 at the top) and axis
% tight.
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
%        Use the per-event-type line-plot format only when the number of event
%        types of interest is <= k (default 12).
%
% :See also:
%   - diagnostics, fmri_glm_design_matrix.plot, drawbox
%
% ..
%    2026 - Renders the object's own basis-convolved regressors per event type
%    (one line per basis function) instead of re-convolving with a canonical
%    HRF, so the actual basis set is respected.
% ..

X = obj.X;
if isempty(X)
    error('glm_map:NoDesign', 'No design matrix available (obj.X is empty).');
end

max_event_types = 12;
wh = find(strcmpi(varargin, 'max_event_types'));
if ~isempty(wh), max_event_types = varargin{wh(1) + 1}; end

% Per-event-type info for session 1 (onsets/durations/names + the X columns
% that hold each condition's basis-function regressors).
ev = local_event_info(obj);
do_panels = ~isempty(ev) && numel(ev) <= max_event_types;

if do_panels
    nconds = numel(ev);
    colors = local_colors(nconds);

    create_figure('glm_map design');

    % Left column: one panel per event type. Right column: full-design heat map.
    for i = 1:nconds
        subplot(nconds, 2, 2 * i - 1);
        local_event_panel(ev(i), colors{i}, i == nconds);
    end

    subplot(nconds, 2, 2:2:2 * nconds);
    local_heatmap(obj, X);
else
    create_figure('glm_map design');
    local_heatmap(obj, X);
end

end % plot_design


% =========================================================================
% Local helpers
% =========================================================================
function ev = local_event_info(obj)
% Build a struct array (one per event type of interest, session 1) with:
%   .name, .time (sec), .Xcols ([nscan1 x nbf] regressor block), .ons, .dur
% Returns [] if this is not a built event design.
ev = [];

if isempty(obj.design) || ~isa(obj.design, 'fmri_glm_design_matrix'), return, end
if isempty(obj.design.Sess) || ~isfield(obj.design.Sess(1), 'U') || isempty(obj.design.Sess(1).U), return, end

X = obj.X;
if isempty(X), return, end

U  = obj.design.Sess(1).U;
TR = obj.TR;

ns1 = obj.design.nscan(1);
if isempty(ns1) || ns1 < 1 || ns1 > size(X, 1), ns1 = size(X, 1); end
time = (0:ns1 - 1)' * TR;

% onset/duration units -> seconds
to_sec = 1;
if ~isempty(obj.design.xBF) && isfield(obj.design.xBF(1), 'UNITS') ...
        && any(strcmpi(obj.design.xBF(1).UNITS, {'scans', 'tr', 'trs'}))
    to_sec = TR;
end

% Interest basis-function columns are ordered condition-by-condition at the
% front of X (session 1 block), nbf columns per condition.
col = 0;
nconds = numel(U);
ev = repmat(struct('name', '', 'time', time, 'Xcols', [], 'ons', [], 'dur', []), 1, nconds);

for i = 1:nconds
    nbf = size(obj.design.xBF(min(i, numel(obj.design.xBF))).bf, 2);
    cols = col + (1:nbf);
    col = col + nbf;
    if any(cols > size(X, 2)), cols = cols(cols <= size(X, 2)); end

    nm = U(i).name; if iscell(nm), if isempty(nm), nm = sprintf('Cond%d', i); else, nm = nm{1}; end, end

    o = U(i).ons(:) * to_sec;
    dd = []; if isfield(U(i), 'dur'), dd = U(i).dur; end
    if isempty(dd), dd = 0; end
    dd = dd(:); if isscalar(dd), dd = repmat(dd, numel(o), 1); end
    dd = dd * to_sec;

    ev(i).name  = nm;
    ev(i).Xcols = X(1:ns1, cols);
    ev(i).ons   = o;
    ev(i).dur   = dd;
end
end


function local_event_panel(e, color, is_last)
% One panel: the event type's actual regressor line(s) (one per basis
% function) plus boxes marking each event and its duration.

Xi = e.Xcols;
nbf = size(Xi, 2);

% Regressor lines: solid color, slightly varying line style per basis function
hold on;
styles = {'-', '--', ':', '-.'};
for j = 1:nbf
    plot(e.time, Xi(:, j), 'Color', color, 'LineStyle', styles{mod(j - 1, numel(styles)) + 1}, 'LineWidth', 1.5);
end

% Event boxes along the bottom of the panel
yl = [min([Xi(:); 0]) max([Xi(:); eps])];
if diff(yl) == 0, yl = yl + [-1 1]; end
boxh = 0.12 * diff(yl);
ybot = yl(1) - 1.6 * boxh;
minw = max(0.4 * (e.time(2) - e.time(1)), 0.001);   % give impulses a visible width
for k = 1:numel(e.ons)
    w = max(e.dur(k), minw);
    local_drawbox(e.ons(k), w, ybot, boxh, color);
end

axis tight;
ylim([ybot - 0.4 * boxh, yl(2) + 0.05 * diff(yl)]);
ylabel(e.name, 'Interpreter', 'none', 'Rotation', 0, 'HorizontalAlignment', 'right');
set(gca, 'YTick', []);
% Time axis is always in seconds (onsets/durations converted from the design's
% units, and the regressor time base is (0:nscan-1)*TR).
if is_last, xlabel('Time (seconds)'); else, set(gca, 'XTickLabel', []); end
if nbf > 1
    legend(arrayfun(@(j) sprintf('BF%d', j), 1:nbf, 'UniformOutput', false), ...
        'Location', 'northeast', 'Box', 'off');
end
hold off;
end


function local_drawbox(t, dur, ystart, yheight, color)
% Filled rectangle for one event (borrowed from drawbox.m).
x = [0 1 1 0] * dur + t;
y = [0 0 1 1] * yheight + ystart;
fill(x, y, color, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end


function c = local_colors(n)
% n distinct colors as a cell array.
base = get(groot, 'DefaultAxesColorOrder');
if isempty(base), base = lines(7); end
while size(base, 1) < n, base = [base; base]; end %#ok<AGROW>
c = mat2cell(base(1:n, :), ones(n, 1), 3)';
end


function local_heatmap(obj, X)
% Heat map of the full design matrix (incl. nuisance covariates / intercept).
% Each regressor (column) is scaled to unit L2 norm for display, so that
% regressors with very different magnitudes (e.g. tiny event regressors and a
% constant intercept) are shown on a comparable scale instead of the
% large-norm columns washing the others out to black.
nrm = vecnorm(X, 2, 1);
nrm(nrm == 0) = 1;
Xdisp = X ./ nrm;

imagesc(Xdisp);
colormap(gray);
colorbar;
set(gca, 'YDir', 'reverse');               % observation 0 at the top
axis tight;
xlabel('Regressor');
ylabel('Image / observation');
title(sprintf('Full design matrix, columns scaled to unit norm (%d interest, %d nuisance, %d intercept)', ...
    sum(obj.wh_interest), sum(obj.wh_nuisance), sum(obj.wh_intercept)));

rn = obj.regressor_names;
if ~isempty(rn) && numel(rn) == size(X, 2)
    labels = rn(:)';
    for j = 1:numel(labels)
        if obj.wh_intercept(j),    labels{j} = [labels{j} ' (icpt)'];
        elseif obj.wh_nuisance(j), labels{j} = [labels{j} ' (nuis)'];
        end
    end
    set(gca, 'XTick', 1:size(X, 2), 'XTickLabel', labels, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
end
end
