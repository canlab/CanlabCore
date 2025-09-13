function [spm_imgs] = expand_4d_filenames(imgs, nvols)
% Expand a list of image names to 4-D name format as used by SPM2/5
% Also compatible with SPM8.
%
% :Usage:
% ::
%
%     [spm_imgs] = expand_4d_filenames(imgs, [nvols])
%
% This function runs with either imgs, nvols, or both variables defined
% If the images exist, spm_imgs is returned as a char matrix with each row
% in the following format: 'path/filename.ext,[N]' where [N] is the number
% of the given image.  The space after [N] is padded with blanks to ensure
% that all the image names fit properly.
%
% :Example:
% ::
%
%    spm_imgs = expand_4d_filenames('run03.img', 10)


    spm_imgs = {};
    if(~exist('imgs', 'var') || isempty(imgs)), error('imgs variable is empty in expand_4d_filenames'); end
    imgs = cellstr(imgs);
    num_files = length(imgs);

    if(exist('nvols', 'var'))
        if(length(nvols) == 1)
            nvols = 1:nvols;
        end
        nvols = repmat({nvols}, num_files, 1);
    end

    for i = 1:num_files
        img_name = imgs{i};
        
        if(any(img_name == ','))
            spm_imgs{end+1} = img_name;
            
        elseif(~exist(img_name, 'file'))
            
            spm_imgs(i) = cell(1, 1);

            if ~exist('nvols', 'var')
                % not existing image AND no vol information
                spm_imgs{i}{1} = [];
                
            else
                % no existing images, BUT we have nvols
                % create a list of to-be-created images
                fprintf('No existing file (yet): %s\n', img_name);
                disp('Returning names of to-be-created volumes.')

                for j = nvols{i}
                    curr_frame = [img_name ',' int2str(j)];
                    spm_imgs{i}{end+1} = curr_frame;
                end

            end

            if iscell(spm_imgs{i})
                spm_imgs{i} = strvcat(spm_imgs{i}{:});
            end
            
        else    %existing image, but unknown volume number...
                % .. basically just uses spm_vol
            if(nargin < 2)
                nvols{i} = 1:num_frames(img_name);
            end

            for j = nvols{i}
                curr_frame = [img_name ',' int2str(j)];
                spm_imgs{end+1} = curr_frame;
            end
        end
    end
    
    spm_imgs = char(spm_imgs);
end



% .img files can contain multiple volumes (indexed by .n in spm_vol structure).
% number of images n stored in this file
% % -------------------------------------------------------------------
function [n] = num_frames(img_name)
    
    % This used FSL if available, but I commented it out.
    % %     global FSLDIR;
    % %
    % %     if(isempty(FSLDIR))
    % %         scn_setup();
    % %     end
    % %
    % %     if(exist(sprintf('%s/bin/fslnvols', FSLDIR), 'file'))
    % %         cmd = sprintf('export FSLDIR=%s && . %s/etc/fslconf/fsl.sh && %s/bin/fslnvols "%s" ', FSLDIR, FSLDIR, FSLDIR, img_name);
    % %         [error_status, result] = system(cmd);
    % %         if(error_status)
    % %             error(result);
    % %         else
    % %             n = str2double(result);
    % %         end
    % %     else
    
    switch(spm('ver'))
        case 'SPM2'
            V = spm_vol(img_name);
            
            fp   = fopen(img_name);
            fseek(fp,0,'eof');
            Len  = ftell(fp);
            fclose(fp);
            n    = Len/(prod(V.dim(1:3))*spm_type(V.dim(4),'bits')/8);

        case {'SPM5', 'SPM8','SPM12', 'SPM25'}
            V = spm_vol(img_name);
            n = length(V);

        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end
    %     end
    
end
