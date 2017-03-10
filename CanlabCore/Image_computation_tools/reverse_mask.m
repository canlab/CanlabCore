function reverse_mask(inname,outname)
% Changes 1 and greater's to 0's in an .img file, and vice versa
% NaNs are still NaNs.
%
% :Usage:
% ::
%
%     reverse_mask(inname,outname)

Q = spm_imcalc_ui(inname,outname,'abs(i1) < eps');

return

