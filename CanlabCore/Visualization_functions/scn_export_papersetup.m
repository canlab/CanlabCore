function scn_export_papersetup(minsize)
% :Usage:
% ::
%
%    scn_export_papersetup([opt: min size in pixels, default = 400])
%
% set paper size for current figure so that print to png or tiff looks as
% it should (as it does on-screen)
%
% ..
%    tor wager, aug. 06
% ..

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
