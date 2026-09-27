function out = count_regions(obj, varargin)
% count_regions  Subject-consistency count-table / count-map ACROSS brain regions.
%
% Given an ARRAY of per-region rsm objects (one RSM per parcel -- e.g. the
% .per_parcel_rsms emitted by fmri_data.rsa_parcelwise), runs count_map in each
% region and assembles, for a chosen contrast/cell, how many subjects show the
% effect in that region. Optionally projects the per-region counts back onto
% the brain as a statistic_image so the standard .threshold / .montage / .table
% chain just works -- the region-level analogue of count_map, and the map you
% report as "the effect was present in >=8/11 subjects in these regions".
%
% Multiple-comparison scope is ACROSS REGIONS (as in rsa_parcelwise): the
% group p-values from each region are FDR-corrected across regions, per contrast.
%
% Usage
% -----
%   res = rsa_parcelwise(dat, 'atlas',atl, 'group_by',{'condition','bodysite'}, ...
%                             'subject_var','sub', 'contrasts',spec);
%   out = count_regions(res.per_parcel_rsms, 'Atlas', atl, ...
%             'Contrasts', { struct('within','hot','between','warm','name','hotVsWarm') });
%   out.region_table                                   % long table, one row per region x contrast
%   montage(threshold(out.maps.hotVsWarm, 0.05, 'unc'));   % paint the count-map
%
% Inputs
% ------
%   obj   [nRegions x 1] array of rsm objects (each k x k x N, same conditions).
%
% Optional name-value
% -------------------
%   'Atlas'       atlas object matching obj (num_regions == numel(obj)); enables
%                 the brain count-map(s). Default [] (table only).
%   'RegionLabels' cellstr region names. Default atlas.labels, else 'region_i'.
%   'Contrasts'   cell of contrast specs (as in count_map: {a,b} pairs or
%                 struct('within',a,'between',b,'name',nm)). Each yields one
%                 value per region. Default: overall within-RSM similarity.
%   'Cell'        alternative to Contrasts for a single block/condition cell:
%                 {name_a,name_b} or [i j]. Ignored if Contrasts given.
%   'MapValue'    what to paint into the brain map(s): 'count' (default) |
%                 'proportion' | 'group_t' | 'group_mean'.
%   'Criterion'/'Threshold'/'Tail'/'Alpha'/'Nperm'/'Transform'/'Reduction'
%                 forwarded to count_map (per-subject criterion). Defaults match.
%   'GroupNperm'  permutations for each region's group stats. Default 5000.
%   'doplot'      logical (default false) -- bar of counts across regions.
%
% Output struct `out`
% -------------------
%   .region_table  long table: Region, Contrast, count, n, proportion,
%                  group_mean, group_t, group_p, group_p_corr (across regions), sig.
%   .count_by_region [nRegions x nContrasts] count matrix.
%   .contrast_names  {1 x nContrasts}.
%   .region_labels   {nRegions x 1}.
%   .maps          struct maps.<contrast> = statistic_image (if Atlas given).
%   .figure        handle if doplot.
%
% See also: fmri_data/rsa_parcelwise, rsm/count_map, assign_vals_to_atlas.

nR = builtin('numel', obj);
if nR < 1 || ~isa(obj, 'rsm'), error('rsm:count_regions:input', 'obj must be an array of rsm objects.'); end

p = inputParser;
p.addParameter('Atlas', [], @(x) isempty(x) || isa(x, 'atlas'));
p.addParameter('RegionLabels', {}, @(x) iscell(x) || isstring(x) || isempty(x));
p.addParameter('Contrasts', {}, @iscell);
p.addParameter('Cell', {}, @(x) iscell(x) || isnumeric(x));
p.addParameter('MapValue', 'count', @(x) ischar(x) || isstring(x));
p.addParameter('Criterion', 'sign', @(x) ischar(x) || isstring(x));
p.addParameter('Threshold', 0, @isscalar);
p.addParameter('Tail', 'right', @(x) ischar(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('Nperm', 1000, @(x) isscalar(x) && x >= 100);
p.addParameter('Transform', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('Reduction', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('GroupNperm', 5000, @(x) isscalar(x) && x >= 100);
p.addParameter('doplot', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
o = p.Results;

k = size(obj(1).dat, 1);
contrasts = local_build_contrasts(o.Contrasts, o.Cell, k);
nC = numel(contrasts);
cnames = local_contrast_names(contrasts, obj(1));

fwd = {'Criterion', o.Criterion, 'Threshold', o.Threshold, 'Tail', o.Tail, ...
    'Alpha', o.Alpha, 'Nperm', o.Nperm, 'Transform', o.Transform, ...
    'Reduction', o.Reduction, 'GroupNperm', o.GroupNperm, 'Correction', 'none'};

count = nan(nR, nC); prop = nan(nR, nC); gmean = nan(nR, nC);
gt = nan(nR, nC); gp = nan(nR, nC); Nsub = nan(nR, 1);
for r = 1:nR
    cm = count_map(obj(r), 'Granularity', 'contrasts', 'Contrasts', contrasts, fwd{:});
    T = cm.table; Nsub(r) = cm.n;
    count(r, :) = T.count(1:nC)';  prop(r, :) = T.proportion(1:nC)';
    gmean(r, :) = T.group_mean(1:nC)'; gt(r, :) = T.group_t(1:nC)'; gp(r, :) = T.group_p(1:nC)';
end
N = max(Nsub);

% ---- FDR across regions, per contrast ----------------------------------
gp_corr = nan(nR, nC);
for c = 1:nC, gp_corr(:, c) = local_fdr_bh(gp(:, c)); end
sig = gp_corr < 0.05;

rlab = local_region_labels(o.RegionLabels, o.Atlas, obj, nR);

% ---- long region table -------------------------------------------------
Reg = {}; Con = {}; cnt = []; pr = []; gm = []; tt = []; pp = []; ppc = []; sg = [];
for c = 1:nC
    Reg = [Reg; rlab(:)]; %#ok<AGROW>
    Con = [Con; repmat(cnames(c), nR, 1)]; %#ok<AGROW>
    cnt = [cnt; count(:, c)]; pr = [pr; prop(:, c)]; gm = [gm; gmean(:, c)]; %#ok<AGROW>
    tt = [tt; gt(:, c)]; pp = [pp; gp(:, c)]; ppc = [ppc; gp_corr(:, c)]; sg = [sg; sig(:, c)]; %#ok<AGROW>
end
region_table = table(string(Reg), string(Con), cnt, repmat(N, numel(cnt), 1), pr, gm, tt, pp, ppc, sg, ...
    'VariableNames', {'Region', 'Contrast', 'count', 'n', 'proportion', ...
    'group_mean', 'group_t', 'group_p', 'group_p_corr', 'sig'});

out = struct('region_table', region_table, 'count_by_region', count, ...
    'contrast_names', {cnames}, 'region_labels', {rlab(:)}, 'n', N, 'maps', struct());

% ---- brain count-maps --------------------------------------------------
if ~isempty(o.Atlas)
    valmat = local_mapval(o.MapValue, count, prop, gt, gmean);
    for c = 1:nC
        try
            map = assign_vals_to_atlas(o.Atlas, [], valmat(:, c), 'p_vals', gp(:, c), ...
                'output_type', 'statistic_image', ...
                'dat_descrip', sprintf('count_regions %s = %s', o.MapValue, cnames{c}));
            out.maps.(matlab.lang.makeValidName(cnames{c})) = map;
        catch ME
            warning('rsm:count_regions:paint', 'Could not paint map for %s: %s', cnames{c}, ME.message);
        end
    end
end

if logical(o.doplot), out.figure = local_plot_regions(out, count, cnames, N); end
end


% =========================================================================
function contrasts = local_build_contrasts(cin, cell_sel, k)
if ~isempty(cin), contrasts = cin; return; end
if ~isempty(cell_sel)
    if isnumeric(cell_sel) && numel(cell_sel) == 2
        contrasts = {{cell_sel(1), cell_sel(2)}};
    else
        contrasts = {cell_sel};
    end
    return
end
contrasts = {{1:k}};                                   % default: overall within-RSM similarity
end


function names = local_contrast_names(contrasts, r0)
names = cell(1, numel(contrasts));
for c = 1:numel(contrasts)
    e = contrasts{c};
    if isstruct(e) && isfield(e, 'name'), names{c} = char(e.name);
    elseif isstruct(e) && isfield(e, 'within'), names{c} = sprintf('%s_vs_btw', local_gname(e.within));
    elseif iscell(e) && numel(e) == 2, names{c} = sprintf('%s_vs_%s', local_gname(e{1}), local_gname(e{2}));
    else
        g = e; if iscell(g), g = g{1}; end
        if isnumeric(g) && numel(g) == numel(r0.labels), names{c} = 'overall'; else, names{c} = local_gname(g); end
    end
    if isempty(names{c}), names{c} = sprintf('contrast_%d', c); end
end
end


function s = local_gname(g)
if ischar(g) || isstring(g), s = char(g); else, s = 'idx'; end
end


function rlab = local_region_labels(rin, atl, obj, nR)
if ~isempty(rin), rlab = cellstr(string(rin)); rlab = rlab(:); return; end
if ~isempty(atl) && numel(atl.labels) == nR, rlab = atl.labels(:); return; end
rlab = cell(nR, 1);
for r = 1:nR
    if ~isempty(obj(r).source) && ~strcmpi(obj(r).source, 'custom'), rlab{r} = obj(r).source;
    else, rlab{r} = sprintf('region_%d', r); end
end
end


function V = local_mapval(which, count, prop, gt, gmean)
switch lower(char(which))
    case 'proportion', V = prop;
    case 'group_t',    V = gt;
    case 'group_mean', V = gmean;
    otherwise,         V = count;
end
end


function q = local_fdr_bh(pv)
pv = pv(:); m = sum(isfinite(pv));
q = nan(size(pv));
idx = find(isfinite(pv));
[ps, ord] = sort(pv(idx));
qq = ps .* m ./ (1:numel(ps))';
qq = min(1, flipud(cummin(flipud(qq))));
tmp = nan(numel(idx), 1); tmp(ord) = qq;
q(idx) = tmp;
end


function h = local_plot_regions(out, count, cnames, N)
h = figure('Color', 'w', 'Name', 'rsm count_regions');
bar(count); ylim([0 N]); ylabel('# subjects'); box off;
set(gca, 'XTick', 1:numel(out.region_labels), 'XTickLabel', out.region_labels, 'XTickLabelRotation', 45);
legend(cnames, 'Location', 'best', 'Interpreter', 'none'); legend boxoff;
title(sprintf('subjects with effect per region (n=%d)', N), 'Interpreter', 'none');
end
