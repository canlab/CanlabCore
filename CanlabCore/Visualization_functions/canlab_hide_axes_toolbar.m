function canlab_hide_axes_toolbar(h)
% Hide the per-axes interaction toolbar and disable default interactions.
%
% Removes the "..." overflow / hover toolbar (pan/zoom/rotate/datatips) that
% recent MATLAB releases (e.g. R2026a) show in the corner of every axes. Static
% brain slices and rendered surfaces don't need those controls, and the button
% clutters every panel of a montage / surface figure.
%
% :Usage:
% ::
%
%     canlab_hide_axes_toolbar(ax)      % one or more axes handles
%     canlab_hide_axes_toolbar(fig)     % every axes in a figure
%     canlab_hide_axes_toolbar(gcf)
%
% Accepts axes handles, figure handles, or a mix; figures are expanded to all of
% their axes. Wrapped in try/catch so it is a safe no-op on older MATLAB.
%
% :See also: fmridisplay, montage, surface, canlab_results_fmridisplay

if nargin < 1 || isempty(h), return, end

axs = gobjects(0);
for g = reshape(h, 1, [])
    if ~isgraphics(g), continue, end
    switch get(g, 'type')
        case 'figure'
            axs = [axs; findall(g, 'Type', 'axes')]; %#ok<AGROW>
        case 'axes'
            axs = [axs; g]; %#ok<AGROW>
    end
end

for a = reshape(axs, 1, [])
    try, a.Toolbar.Visible = 'off'; catch, end
    try, disableDefaultInteractivity(a); catch, end
    try, a.Interactions = []; catch, end
end
end
