function mask = mask_create_from_image_set(imgs, outname, atleastN, varargin)
% Take a set of images and create a mask of voxels in which at least N
% subjects have valid (not exactly zero, non NaN) data.
%
% :Usage:
% ::
%
%     mask = mask_create_from_image_set(imgs, outname, atleastN, ['sum'])
%
% This makes a useful results mask for a set of images, i.e., in a
% group analysis.
%
% :Optional: 'sum' input writes the sum image instead of the mask image,
% so that the values in the image reflect the number of input images
% with valid values.
%
% compatible with SPM5 and above only!
%
% :Examples:
% ::
%
%    mask_create_from_image_set(EXPT.SNPM.P{1}, 'mask_all_valid.img');
%
%    imgs = filenames('vascular_mask_*img');
%    mask_create_from_image_set(imgs, 'vascular_group_sum.img', 6, 'sum');
%
% ..
%    Tor Wager, April 2, 2008
%    Edit Sept 13, 2008: Sum image output option
% ..

    write_sum = 0;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case {'sum', 'sumimage'}, write_sum = 1; 
                
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    imgs = strvcat(imgs); % enforce string
    
    if nargin < 3 || isempty(atleastN), atleastN = size(imgs, 1); end
    
    spm_defaults
    
    disp('Mapping images.');
    V = spm_vol(imgs);
    
    disp('Checking dimensions and spaces.');
    anybad = iimg_check_volinfo(V(1), V);
    if anybad, error('Exiting.'); end
    
    disp('Reading images.');
    dat = spm_read_vols(V);
    
    % number of subjects for which each vox is valid  
    sumimg = sum(abs(dat) > eps & ~isnan(dat), 4);
    
    % threshold
    mask = double(sumimg >= atleastN);
    
    
    % write output
    if ~isempty(outname)
        
        
        outV = struct( 'fname', outname, 'mat', V(1).mat, 'dim', V(1).dim, 'dt', V(1).dt, 'pinfo', [1 0 0]', 'n', [1 V(1).n(2)]  );
        outV.descrip = sprintf('At least %3.0f out of %3.0f images valid', atleastN, size(imgs, 1));
        
        fprintf('Writing: %s\n', outV.fname);
        
        if write_sum
            spm_write_vol(outV, sumimg);
        else
            spm_write_vol(outV, mask);
        end
    end
    
end
