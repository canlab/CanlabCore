function voldata = iimg_reconstruct_3dvol(dat, volInfo, varargin)
% Reconstruct a 3-D volume from a dat index list
%
% :Usage:
% ::
%
%     voldata = iimg_reconstruct_3dvol(dat, volInfo, [optional args])
%
% :Optional Inputs: if entered, will write .img file to disk
%
%   **'outname':**
%        followed by output name
%
%   **'descrip':**
%        followed by description for img file
%
%   **'slice':**
%        followed by slice number of single-slice data in image
% 
% THIS FUNCTION IS DEPRECATED. USE IIMG_RECONSTRUCT_VOLS.M'
% WHICH CAN DEAL WITH 4-D VOLUMES AS WELL


disp('THIS FUNCTION IS DEPRECATED. USE IIMG_RECONSTRUCT_VOLS.M');

% Do not force spm scaling factors and get rid of private, which may
% contain information from other images/irrelevant info for this...
if isfield(volInfo, 'pinfo'), volInfo = rmfield(volInfo, 'pinfo'); end
if isfield(volInfo, 'private'), volInfo = rmfield(volInfo, 'private'); end

    is_writing_file = 0;
    is_single_slice = 0;
    slice_number = NaN;
    descrip = 'Created by iimg_reconstruct_3dvol';

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'outname'
                    outname = varargin{i+1};
                    is_writing_file = 1;

                case 'descrip'
                    descrip = varargin{i+1};

                case 'slice'
                    is_single_slice = 1;
                    slice_number = varargin{i+1};
            end
        end
    end

    % check data and make sure it conforms to full-length index vector
    % adjust if only in-mask voxels are in data
    dat = check_dat(dat, volInfo, is_single_slice, slice_number);

    if is_single_slice
        voldata = reshape(dat, volInfo.dim(1:2));
    else
        voldata = reshape(dat, volInfo.dim(1:3));
    end

    if is_writing_file
        if is_single_slice
          write_file(voldata, volInfo, outname, descrip, slice_number);  
        else
        write_file(voldata, volInfo, outname, descrip);
        end
    end
    
end  % End Main Function


% ------------------------------------------------------------------
% * Check and make sure data is in full list; adjust if mask-only
% ------------------------------------------------------------------
function dat = check_dat(dat, volInfo, is_single_slice, varargin)
    % last arg is slice no., if applicable

    if ~isfield(volInfo, 'nvox')
        error('volInfo.nvox and other fields not found.  Use iimg_read_image.m to prepare volInfo.');
    end

    ndat = length(dat);

    if is_single_slice
        % data should either equal prod(volInfo(1:2)) or fewer elements, in
        % which case one slice will be reconstructed
        nslice = prod(volInfo.dim(1:2));
        z = volInfo.xyzlist(:,3);
        slice_number = varargin{1};

        if ndat ~= nslice
            % get indices in full volume of this slice
            st = 1 + (slice_number - 1) .* nslice;  % starting index of this slice
            en = slice_number .* nslice;

            tmpdat = zeros(nslice, 1); % full volume, so we know how voxels map into 3-D array
            tmpdat(volInfo.image_indx(st:en)) = dat;
            dat = tmpdat;
        end

        % if data already includes all voxels, do nothing, otherwise assume data is only in-mask, and reconstruct
    elseif ndat ~= volInfo.nvox
        if ~isfield(volInfo, 'n_inmask')
            error('volInfo.n_inmask not found.  Use iimg_read_img.m with extended output flag to prepare volInfo.');
        end

        if ndat ~= volInfo.n_inmask
            error('volume info in struct does not seem to match number of voxels in index image.');
        else
            tmpdat = zeros(volInfo.nvox, 1);
            tmpdat(volInfo.image_indx) = dat;
            dat = tmpdat;
        end
    end
end

% ------------------------------------------------------------------
% Write full image or slice, if last arg. entered
% ------------------------------------------------------------------
function write_file(voldata, volInfo, outname, descrip, varargin)
    [volInfo.fname, volInfo.n(1)] = make_img_filename(outname);
    volInfo.descrip = descrip;

    if length(varargin) > 0
        slice_number = varargin{1};
        % SPM5 uses .private field in spm_write_plane...make sure we have
        % correct .private loaded from file by re-mapping volume.
        spm_write_plane(spm_vol(volInfo.fname), voldata, slice_number);
    else
        warning('off'); % empty images return many warnings
        spm_write_vol(volInfo, voldata);
        warning('on');
    end
end

function [name, n] = make_img_filename(name)
[d, name, ext] = fileparts(name);
name(name == '.') = '_';

% changed nov 07 to preserve commas for multi-volume images.
% handle the comma to index volume: return as n and take it off the
% filename
t = find(ext==',');

n = 1;
if ~isempty(t)
    t = t(1);
    n1 = ext((t+1):end);
    if ~isempty(n1),
        n   = str2num(n1);
        ext = ext(1:(t-1));
    end
end

name = fullfile(d, [name ext]);

if ~strcmp(name(end-3:end), '.img') && ~strcmp(name(end-3:end), '.nii')
    name = [name '.img'];
end
end
