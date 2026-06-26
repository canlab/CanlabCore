function han = plot(obj, varargin)
% plot Quick QC plots for an fmri_surface_data object.
%
% :Usage:
% ::
%     han = plot(obj)
%
% Surface analogue of fmri_data.plot: a compact QC panel with (1) a histogram of
% grayordinate values, (2) the per-map mean/standard-deviation across
% grayordinates, and (3) a native surface render of the mean map. Geometry-
% agnostic summaries operate directly on .dat.
%
% :Optional Inputs:
%   **'norender':** skip the surface render panel (faster; no mesh load).
%
% :Outputs:
%   **han:** struct with handles (.figure, and .surface for the render).
%
% :See also: surface, descriptives, fmri_surface_data

dorender = ~any(strcmpi(varargin, 'norender'));

d = double(obj.dat);
nGray = size(d, 1);
nMaps = size(d, 2);

fig = figure('Color', 'w', 'Name', 'fmri_surface_data QC');

% (1) Histogram of all values
subplot(2, 2, 1);
vals = d(:); vals = vals(~isnan(vals));
histogram(vals, 100);
title(sprintf('Value histogram (%d grayordinates x %d maps)', nGray, nMaps));
xlabel('value'); ylabel('count'); axis tight;

% (2) Per-map mean +/- sd
subplot(2, 2, 2);
mu = mean(d, 1, 'omitnan');
sd = std(d, 0, 1, 'omitnan');
errorbar(1:nMaps, mu, sd, 'o-', 'LineWidth', 1);
title('Per-map mean \pm sd across grayordinates');
xlabel('map'); ylabel('value'); xlim([0.5 nMaps + 0.5]); grid on;

% (3) Fraction nonzero per map (coverage)
subplot(2, 2, 3);
cov = mean(d ~= 0 & ~isnan(d), 1);
bar(1:nMaps, cov);
title('Fraction nonzero per map'); xlabel('map'); ylabel('fraction'); ylim([0 1]);

han = struct('figure', fig, 'surface', []);

% (4) Surface render of the mean map
if dorender && ~isempty(obj.brain_model) ...
        && any(cellfun(@(m) strcmp(m.type,'surf'), obj.brain_model.models))
    try
        mobj = mean(obj);
        ax = subplot(2, 2, 4);
        geom = load_surface_geom(obj.surface_space, 'inflated');
        hp = patch('Parent', ax, 'Faces', geom.faces_lh, 'Vertices', geom.vertices_lh, ...
            'EdgeColor', 'none', 'FaceColor', [.5 .5 .5], 'Tag', 'left');
        axis(ax, 'off', 'image', 'vis3d'); view(ax, 270, 0);
        try, lightRestoreSingle(ax); catch, camlight(ax); end %#ok<NOCOM>
        material(ax, 'dull');
        render_on_surface(mobj, hp);
        title(ax, 'Mean map (left lateral)');
        han.surface = hp;
    catch err
        warning('fmri_surface_data:plot:render', 'Surface render skipped: %s', err.message);
    end
end
end
