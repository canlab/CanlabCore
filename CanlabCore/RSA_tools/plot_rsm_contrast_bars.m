function h = plot_rsm_contrast_bars(data, varargin)
% plot_rsm_contrast_bars  Bar plot of per-subject RSA contrast values.
%
% Wraps barplot_columns + plotWithinGroup with the styling conventions used
% across the Sun et al. workflow notebooks (Figs 7D-F): bar plus per-subject
% lines, optional pairwise-significance markers, within-subject paired t-test
% defaults.
%
% Usage
% -----
%   % From an rsm: pass a cells_table or contrast spec
%   T = R.cells_table({'hot','warm','imagine'});
%   h = plot_rsm_contrast_bars(T);
%
%   % Or directly from a per-subject matrix or cell array
%   h = plot_rsm_contrast_bars(matrix, 'names', {'Hot','Warm','Imagine'});
%   h = plot_rsm_contrast_bars({v1, v2, v3}, 'names', {'Hot','Warm','Imagine'});
%
% Optional name-value
% -------------------
%   'names'         cellstr of bar labels (default: from table columns or auto)
%   'colors'        cell of RGB triplets, one per bar (default: hsv)
%   'within_lines'  true (default) -- draw per-subject connecting lines via
%                   plotWithinGroup
%   'within_colors' [n_subjects x 3] RGB matrix for the within-subject lines
%   'pairwise'      true (default) -- show pairwise significance bars
%   'ylabel'        default 'Fisher r-to-z'
%   'title'         default ''
%   'ylim'          default [] (auto)
%   'face_alpha'    default 1
%
% Output
% ------
%   h  struct with handles from barplot_columns

p = inputParser;
p.addParameter('names',         {},     @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('colors',        {},     @(x) iscell(x) || isnumeric(x) || isempty(x));
p.addParameter('within_lines',  true,   @(x) islogical(x) || isnumeric(x));
p.addParameter('within_colors', [],     @isnumeric);
p.addParameter('pairwise',      true,   @(x) islogical(x) || isnumeric(x));
p.addParameter('ylabel',        'Fisher r-to-z', @(x) ischar(x) || isstring(x));
p.addParameter('title',         '',     @(x) ischar(x) || isstring(x));
p.addParameter('ylim',          [],     @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
p.addParameter('face_alpha',    1,      @isnumeric);
p.parse(varargin{:});
opt = p.Results;

% Coerce data to a cell of column vectors and pull names
if istable(data)
    cols = data.Properties.VariableNames;
    cell_data = cellfun(@(c) data.(c), cols, 'UniformOutput', false);
    if isempty(opt.names), opt.names = cols; end
elseif ismatrix(data) && ~iscell(data)
    cell_data = num2cell(data, 1);
    if isempty(opt.names)
        opt.names = arrayfun(@(i) sprintf('condition_%d', i), 1:size(data,2), 'UniformOutput', false);
    end
elseif iscell(data)
    cell_data = data;
    if isempty(opt.names)
        opt.names = arrayfun(@(i) sprintf('condition_%d', i), 1:numel(cell_data), 'UniformOutput', false);
    end
else
    error('plot_rsm_contrast_bars:badData', ...
        'data must be a table, matrix, or cell of vectors.');
end

names = cellstr(opt.names);
n_bars = numel(cell_data);

% Default colors
if isempty(opt.colors)
    cmap = hsv(n_bars);
    colors = num2cell(cmap, 2);
else
    colors = opt.colors;
end
if isnumeric(colors), colors = num2cell(colors, 2); end

% Build barplot_columns args
args = {'names', names, 'color', colors, 'noviolin', 'noind', 'nofig'};
if opt.pairwise, args = [args, {'pairwisetest'}]; end

h = barplot_columns(cell_data, args{:});
hold on;

% Style bars
for k = 1:numel(h.bar_han)
    h.bar_han{k}.LineWidth = 3;
    h.bar_han{k}.FaceAlpha = opt.face_alpha;
end

% Within-subject lines
if opt.within_lines && exist('plotWithinGroup', 'file') == 2
    wg_args = {'PlotHalfViolin', false, 'PlotBox', false};
    if ~isempty(opt.within_colors)
        wg_args = [wg_args, {'IndColors', opt.within_colors}];
    else
        % Default: hsv per subject
        n_subj = numel(cell_data{1});
        wg_args = [wg_args, {'IndColors', hsv(n_subj)}];
    end
    try
        plotWithinGroup(cell_data, names, wg_args{:});
    catch
        % silently skip if plotWithinGroup signature differs
    end
end

% Zero line, labels
if exist('hline', 'file') == 2
    try, hline(0, 'k-'); catch, end
end
ylabel(char(opt.ylabel));
xlabel('');
if ~isempty(char(opt.title)), title(char(opt.title), 'Interpreter','none'); end
if ~isempty(opt.ylim), ylim(opt.ylim); end
xtickangle(0);

hold off;

end
