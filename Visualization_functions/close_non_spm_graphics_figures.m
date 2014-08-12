function close_non_spm_graphics_figures
%close_non_spm_graphics_figures
%
% closes all figure windows that are not the SPM orthviews window


fh = findobj('Type','Figure');
spmf = findobj('Tag','Graphics');
if ~isempty(spmf), fh(fh - spmf == 0) = []; end
close(fh)

return