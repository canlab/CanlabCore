function h = image_scatterplot(tmap1, tmap2, varargin)
% IMAGE_SCATTERPLOT Scatter/contour comparison of two CANlab image objects (t-values)
%
% :Usage:
% ::
%     h = image_scatterplot(tmap1, tmap2, ['pvaluebox', p], ['df', df], ...
%                           ['colorpoints'], ['contour'], ...
%                           ['nbins', nbins], ['title', titlestr]);
%
% :Inputs:
%   **tmap1, tmap2**
%       CANlab image objects with a .dat field (e.g., fmri_data/statistic_image).
%       This function plots tmap1.dat (x-axis) against tmap2.dat (y-axis).
%
% :Optional Inputs (name/value pairs and flags):
%   **'pvaluebox', p**
%       Draw a thin black dashed square corresponding to the two-tailed threshold
%       p < p. The threshold is computed as |t| >= tinv(1 - p/2, df).
%       If df is not available, the function falls back to a normal approximation
%       |z| >= norminv(1 - p/2).
%
%   **'df', df**
%       Degrees of freedom used to convert p-values to a t-threshold for 'pvaluebox'.
%       If omitted, this function attempts to infer df from the input objects (common
%       fields include .dfe or .df). If inference fails, a normal approximation is used.
%
%   **'colorpoints'**
%       Color-codes a subset of points:
%         - Dark pink: x > 0 and y > x  OR  x < 0 and y < x
%           (map2 is more extreme than map1 in the same sign as map1)
%         - Blue: y > 0 and x > y  OR  y < 0 and x < y
%           (map1 is more extreme than map2 in the same sign as map2)
%       All other points are plotted in a default neutral color.
%
%   **'contour'**
%       Plot a contour plot of a 2D histogram (density) instead of individual points.
%
%   **'nbins', nbins**
%       Number of bins per dimension for the 2D histogram used in 'contour' mode.
%       Default: 60
%
%   **'title', titlestr**
%       Figure title string. Default: 't-value comparison'
%
% :Outputs:
%   **h**
%       Struct of handles and useful outputs:
%         h.fig, h.ax, h.scatter, h.contour, h.identity, h.pbox, h.corrtext
%         h.r (Pearson correlation), h.n (number of valid points)
%
% :Details:
%   - Axis limits are set dynamically from the data, symmetric around zero, and the
%     same limits are applied to x and y (square axes).
%   - The Pearson correlation between x and y is computed and displayed on the plot
%     in black 18 pt font.
%
% :Examples:
% ::
%     % Basic scatter
%     image_scatterplot(t_obj_temp_2sss, t_obj_temp_fitlme);
%
%     % Add a p-value threshold box (two-tailed) using inferred df (or pass df explicitly)
%     image_scatterplot(t_obj_temp_2sss, t_obj_temp_fitlme, 'pvaluebox', 0.001, 'df', 120);
%
%     % Color points by relative extremity
%     image_scatterplot(t_obj_temp_2sss, t_obj_temp_fitlme, 'colorpoints');
%
%     % Contour (density) instead of points
%     image_scatterplot(t_obj_temp_2sss, t_obj_temp_fitlme, 'contour', 'nbins', 80);
%
% ..
%    Author: ChatGPT (adapted to CANlab-style inputParser + documentation_template.m style)
% ..

% -------------------------------------------------------------------------
% Input parsing (CANlab-style)
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'tmap1', @(x) isstruct(x) || isobject(x));
addRequired(p, 'tmap2', @(x) isstruct(x) || isobject(x));

addParameter(p, 'pvaluebox', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x > 0 && x < 1));
addParameter(p, 'df', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x > 0));
addParameter(p, 'nbins', 60, @(x) isscalar(x) && isnumeric(x) && x >= 10);

addParameter(p, 'title', 't-value comparison', @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(p, 'colorpoints', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'contour', false, @(x) islogical(x) && isscalar(x));

% Support flag-style usage: 'colorpoints', 'contour'
% (inputParser does not natively treat bare strings as flags unless paired)
varargin = local_coerce_flags(varargin, {'colorpoints', 'contour'});

parse(p, tmap1, tmap2, varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Extract data
% -------------------------------------------------------------------------
x = local_get_dat_vector(tmap1);
y = local_get_dat_vector(tmap2);

% Ensure column vectors
x = x(:); y = y(:);

% Drop non-finite pairs
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);

h = struct();
h.n = numel(x);

% -------------------------------------------------------------------------
% Create figure/axes
% -------------------------------------------------------------------------
try
    h.fig = create_figure(opt.title);
catch
    h.fig = figure('Name', opt.title);
end
set(h.fig, 'Color', 'w');
set(h.fig, 'DefaultAxesFontSize', 18);
h.ax = gca;
hold(h.ax, 'on');

% Dynamic symmetric limits shared across x/y
maxAbs = max(abs([x; y]));
if isempty(maxAbs) || maxAbs == 0
    maxAbs = 1;
end
lims = [-maxAbs maxAbs];

% -------------------------------------------------------------------------
% Plot: contour or scatter
% -------------------------------------------------------------------------
h.scatter = [];
h.contour = [];

if opt.contour
    nbins = opt.nbins;

    % Use symmetric edges so the density grid matches the axes
    edges = linspace(lims(1), lims(2), nbins + 1);

    % histcounts2: rows correspond to x-bins, columns to y-bins
    N = histcounts2(x, y, edges, edges);

    % Convert bin edges to centers for contour
    ctrs = (edges(1:end-1) + edges(2:end)) ./ 2;
    [Xc, Yc] = meshgrid(ctrs, ctrs);

    % Note transpose so contour axes align as expected
    h.contour = contour(h.ax, Xc, Yc, N', 'LineWidth', 1);
else
    if opt.colorpoints
        % Base points (neutral)
        neutral = [0.35 0.35 0.35];
        h.scatter = plot(h.ax, x, y, '.', 'Color', neutral);

        % Highlight rules
        dark_pink   = [0.9 0.00 0.6];
        green = [0.20 0.7 0.2];

        idx_pink =  (x > 0 & y > x) | (x < 0 & y < x);
        idx_orng =  (y > 0 & x > y) | (y < 0 & x < y);

        plot(h.ax, x(idx_pink), y(idx_pink), '.', 'Color', dark_pink);
        plot(h.ax, x(idx_orng), y(idx_orng), '.', 'Color', green);
    else
        h.scatter = plot(h.ax, x, y, '.');
    end
end

% Identity line
h.identity = plot(h.ax, lims, lims, '--', 'Color', 'k');

% Labels
xlabel(h.ax, 't-values (map 1)');
ylabel(h.ax, 't-values (map 2)');

% Axis formatting
xlim(h.ax, lims);
ylim(h.ax, lims);
axis(h.ax, 'square');

% -------------------------------------------------------------------------
% p-value threshold box (two-tailed)
% -------------------------------------------------------------------------
h.pbox = [];
if ~isempty(opt.pvaluebox)
    pthr = opt.pvaluebox;
    df = opt.df;
    if isempty(df)
        df = local_infer_df(tmap1, tmap2);
    end

    tthr = local_p_to_threshold(pthr, df);

    % Draw square at +/- tthr
    xx = [-tthr  tthr  tthr -tthr -tthr];
    yy = [-tthr -tthr  tthr  tthr -tthr];
    h.pbox = plot(h.ax, xx, yy, '--', 'Color', 'k', 'LineWidth', 1);

    % -------------------------------------------------------------------------
    % Override colors inside p-value box (if both pvaluebox and colorpoints)
    % -------------------------------------------------------------------------
    if opt.colorpoints && ~isempty(opt.pvaluebox)

        % Points inside the two-tailed p-value box
        inside_box = abs(x) < tthr & abs(y) < tthr;

        % Plot over existing points in gray
        plot(h.ax, x(inside_box), y(inside_box), '.', ...
            'Color', [0.5 0.5 0.5]);
    end

end

% -------------------------------------------------------------------------
% Correlation label
% -------------------------------------------------------------------------
if h.n >= 2
    r = corr(x, y, 'Rows', 'complete');
else
    r = NaN;
end
h.r = r;

% Place text near top-left inside axes
tx = lims(1) + 0.05 * range(lims);
ty = lims(2) - 0.08 * range(lims);
h.corrtext = text(h.ax, tx, ty, sprintf('r = %.3f', r), ...
    'Color', 'k', 'FontSize', 22, 'FontWeight', 'normal', 'VerticalAlignment', 'top');

set(gca, 'FontSize', 22)

% -------------------------------------------------------------------------
% Nested/local helpers
% -------------------------------------------------------------------------
end


% ========================================================================
% Helper functions (local)
% ========================================================================

function v = local_get_dat_vector(obj)
% Robustly extract .dat as a numeric vector
if isobject(obj) || isstruct(obj)
    if isfield(obj, 'dat')
        v = obj.dat;
    else
        try
            v = obj.dat;
        catch
            error('image_scatterplot:NoDatField', 'Input object must contain a .dat field.');
        end
    end
else
    error('image_scatterplot:BadInput', 'Inputs must be CANlab-like objects with a .dat field.');
end

if ~isnumeric(v)
    error('image_scatterplot:DatNotNumeric', '.dat must be numeric.');
end
end


function df = local_infer_df(obj1, obj2)
% Try common CANlab/statistic_image df fields; return [] if unknown
df = [];

    df1 = obj1.dfe;
    df2 = obj2.dfe;

    if ~isempty(df1) && isfinite(df1) && df1 > 0
        df = df1; return;
    end
    if ~isempty(df2) && isfinite(df2) && df2 > 0
        df = df2; return;
    end
end


function thr = local_p_to_threshold(p_two_tailed, df)
% Convert two-tailed p to |t| threshold; if df unknown, use normal approx
if isempty(df)
    % Normal approximation
    % thr = abs(norminv(p_two_tailed / 2, 0, 1));
    thr = abs(norminv(1 - p_two_tailed / 2, 0, 1));
else
    thr = tinv(1 - p_two_tailed / 2, df);
end

thr = abs(thr);

end


function args = local_coerce_flags(args, flagNames)
% Turn bare flags like 'contour' into name/value pairs: 'contour', true
% Only for known flags.
i = 1;
out = {};
while i <= numel(args)
    if ischar(args{i}) || (isstring(args{i}) && isscalar(args{i}))
        key = char(args{i});
        if any(strcmpi(key, flagNames))
            out = [out, {key, true}]; %#ok<AGROW>
            i = i + 1;
            continue
        end
    end
    out = [out, args(i)]; %#ok<AGROW>
    i = i + 1;
end
args = out;
end