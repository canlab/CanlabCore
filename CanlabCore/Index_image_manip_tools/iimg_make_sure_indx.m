function dat = iimg_make_sure_indx(inputarg)
% :Usage:
% ::
%
%     dat = iimg_make_sure_indx(inputarg)
%
% Make sure input is an index list, and convert if not

[m,n] = size(inputarg);

if isstr(inputarg)
    % it's a file name
    [V,dat] = iimg_read_image(inputarg);

elseif n > 1
    % it's a matrix
    dat = inputarg(:);
    
else
    % it's an index vector
    % do nothing
    dat = inputarg;
end

return

    


