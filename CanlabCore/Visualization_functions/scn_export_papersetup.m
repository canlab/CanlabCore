function scn_export_papersetup(minsize)
% scn_export_papersetup Set paper size of current figure for export to image files.
%
% :Usage:
% ::
%
%     scn_export_papersetup([minsize])
%
% Set the paper size for the current figure so that printing to PNG or
% TIFF (e.g., via saveas, print, or export_fig) produces output that
% matches what is shown on-screen. The figure aspect ratio is preserved,
% and the smaller dimension of the figure is scaled to minsize
% pixels (default 400).
%
% ..
%    tor wager, aug. 06
% ..
%
% :Inputs:
%
%   **minsize:**
%        Optional. Minimum size in pixels for the smaller of the figure's
%        width and height. Default = 400. The larger dimension is scaled
%        proportionally so the original aspect ratio is preserved.
%
% :Outputs:
%
%   None. The function modifies the 'PaperUnits', 'PaperPosition', and
%   'PaperType' properties of the current figure (gcf) in place.
%
% :Examples:
% ::
%
%     % Set up paper size and save the current figure as a PNG
%     scn_export_papersetup(400);
%     saveas(gcf, 'my_figure.png');
%
%     % Larger paper size for higher-resolution exports
%     scn_export_papersetup(800);
%     print(gcf, '-dpng', '-r300', 'my_figure_hires.png');
%
% :See also:
%   - saveas
%   - print
%   - canlab_results_fmridisplay

if nargin < 1, minsize = 400; end

% make sure that min size of fig. is 400 pixels
sz = get(gcf,'Position');
sz = sz(3:4);   % width and height
szratio = max(sz) ./ min(sz);   % max to min size
wh = find(sz == min(sz));
sz(wh) = minsize;
wh = find(sz == max(sz));
sz(wh) = minsize * szratio;

% set paper size using current screen
set(gcf,'PaperUnits','points','PaperPosition',[0 0 sz],'PaperType','usletter');
return
