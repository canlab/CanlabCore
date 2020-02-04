function cm = spm_orthviews_change_colormap(lowcolor, hicolor, varargin)
% :Usage:
% ::
%
%    cm = spm_orthviews_change_colormap(lowcolor, hicolor, midcolor1, midcolor2, etc.)
%
% Create a new split colormap of your choosing and apply it to the
% spm_orthviews figure.
%
% :Examples:
% ::
%
%    spm_orthviews_change_colormap([.2 .2 .6], [1 1 0]);  % slate to yellow
%    spm_orthviews_change_colormap([.9 .5 .2], [1 1 0]);  % orange to yellow
%    spm_orthviews_change_colormap([.8 .1 .1], [1 1 0], [.9 .6 .1]);  %red to orange to yellow
%    spm_orthviews_change_colormap([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
%    cm = spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 .5 1], [0 .5 .5], ...
%    [0 1 .5], [0 1 0], [.5 1 0]);
%    cm = spm_orthviews_change_colormap([0 0 1], [1 1 0], [.5 0 1], [.5 .5 1], ... 
%    [1 .5 1], [1 .5 .5], [1 .5 0]);
%
% ..
%    tor wager, sept. 2007
% ..

myfig = findobj('Tag','Graphics'); % colormap
axishandles = findobj(myfig, 'Type', 'Axes');

newcm = colormap_tor(lowcolor, hicolor, varargin{:});
wh = round(linspace(1, size(newcm, 1), 64));
cm = [gray(64) ; newcm(wh, :)];


% specify split color-map of 128 values; first 64 are gray, last 64 have
% colors for activations; this is specified by spm_orthviews

% Don't set directly, because legends will not draw right.
% Use spm_figure. % Redraw color bar

spm_figure_canlab('Colormap',cm)

% for i = 1:length(axishandles)
%     colormap(axishandles(i), cm);
% end


end

% 
% newcm = NaN .* zeros(64, 3);  %initialize
% n = 64;
% 
% % if nargin < 3
% %     % no mid-color
% %     newcm = [linspace(lowcolor(1), hicolor(1), n)' linspace(lowcolor(2), hicolor(2), n)' linspace(lowcolor(3), hicolor(3), n)' ];
% % else
% %     % we have a mid-color
% %     newcm1 = [linspace(lowcolor(1), midcolor(1), n/2)' linspace(lowcolor(2), midcolor(2), n/2)' linspace(lowcolor(3), midcolor(3), n/2)' ];
% %     newcm2 = [linspace(midcolor(1), hicolor(1), n/2)' linspace(midcolor(2), hicolor(2), n/2)' linspace(midcolor(3), hicolor(3), n/2)' ];
% %
% %     newcm = [newcm1; newcm2];
% % end
