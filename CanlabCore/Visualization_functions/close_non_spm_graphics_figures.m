function close_non_spm_graphics_figures()

% :Usage:
% ::
%
%    close_non_spm_graphics_figures()
%
% purpose:  closes all figure windows that are not the SPM orthviews window
%
% :Input:
%
%   
%
% :Output:
%
%   
% 


fh = findobj('Type','Figure');
spmf = findobj('Tag','Graphics');
if ~isempty(spmf), fh(fh - spmf == 0) = []; end
close(fh)

return
