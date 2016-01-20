function iimg_write_images(dat, volInfo, outnames)
% Write a series of Analyze images given inputs.
%
% :Usage:
% ::
%
%     iimg_write_images(dat,volInfo,outnames)
%
% :Inputs:
%
%   **dat:**
%        voxels x images matrix of image data in index format
%
%   **volInfo:**
%        spm-style info structure; see iimg_read_img
%
%   **outnames:**
%        a string matrix (or cell array) of output filenames, or the name of a 4-D file to create 

    
    if isempty(outnames)
        error('You must specify a string array of output names, or a single 4-D image file name');
    end

    if iscell(outnames)
        outnames = char(outnames{:});
    end

    n = size(dat,2);
    if size(outnames,1) ~= n
        error('Number of images (columns of dat) must == number of output img names.');
    end

    if volInfo.nvox ~= size(dat,1)
        error('Dims of dat (vox x images) and volInfo do not match.');
    end

    is4d = 0;
    if size(dat, 2) > 1 && size(outnames, 1) == 1, is4d = 1; end
    
    if is4d
        iimg_reconstruct_vols(dat,volInfo,'outname',deblank(outnames));
        
    else
        for i = 1:n
            %iimg_reconstruct_3dvol(dat(:,i),volInfo,'outname',deblank(outnames(i,:)));
            % This handles volume creation, etc.  It also works for a 4-D
            % image (but all volumes must be created at once)
            iimg_reconstruct_vols(dat(:,i),volInfo,'outname',deblank(outnames(i,:)));
        end
    end
    
end
