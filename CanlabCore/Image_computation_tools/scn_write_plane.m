function V = scn_write_plane(filenames_or_V, dat, wh_slice, varargin)
% Images can be 3D or 4D (4D is possible with SPM2+ using ,xx indexing, or SPM5+).
%
% :Usage:
% ::
%
%     V = scn_write_plane(filenames_or_V, dat, wh_slice, [exampleV])
%
% :Inputs:
%
%   **ofilenames_or_V:**
%        is string matrix or cell array of names, or structures created with spm_vol(names)
%
%   **dat:**
%        is 3-D array of data, vox x vox x slices, with data on 3rd dim being
%        written to separate image files
%
%   **wh_slice:**
%        slice number (voxel space)
%
%   **exampleV:**
%        is optional, but required if filenames are passed in.
%
% Write a plane, given filenames or spm_vol structures and 4-D data with
% the slice to write.
% 
% SPM2/5 compatible, and creates image if necessary from file name and
% example V structure.
%
% ::
%
%    P = char({'one.img', 'two.img', 'three.img'})
%    dat = randn(64, 64, 3);
%    Vout = scn_write_plane(P, dat, wh_slice, V)
%    spm_image('init', Vout(1).fname);
%
%
% Notes on SPM5 and why we need to do what we do the way we do it:
% Write data...in SPM5, by accessing file_array object in
% V.private directly (with spm_write_plane. spm_write_vol 
% first creates NIFTI obj/file array and then uses spm_write_plane.)
%
% The file_array object points to data in
% the actual file, so when values are assigned to the array
% object, they are written directly in the file.
% These values depend on the offset and slope values (spm
% scaling factors!) in V.private.dat as well, so care must be
% taken to assign data to a file_array object with the correct
% scaling factors.  This is why it is better to load the
% structure one wants to write to with spm_vol first, so that
% the name in V.private.dat.fname will be correct, and
% V.private.dat.scl_slope and ...inter will be correct as well.
% Vout = spm_vol(V(i));  % loads correct .private info from V.fname
%
% in SPM5, this simply assigns data in Vout.private.dat
% in SPM2, it does something different, but should be
% compatible, since spm_vol was used above...
% ::
%
%    spm_write_plane(Vout, fsl(:, :, i), slicei);
%    if isfield(V, 'dt'), Vout.dt = V.dt; end      % SPM5 only
%    if isfield(V, 'n'), Vout.n = V.n;  end        % SPM5 only
%        Vout = spm_create_vol(Vout);
%    else
%        Vout = spm_vol(V(i).fname);
%    end
%    spm_write_plane(Vout, fsl(:, :, i), slicei);
%
% ..
%    Tor Wager, March 2008
% ..


% ..
%    Make sure all images exist, and create if not
%    Return V structure of spm_vol mapped image structures
% ..

if ~isstruct(filenames_or_V)
    % Assume we have a string matrix of file names
    % Make sure they exist...
    % create volumes, if necessary

    
    % Make sure char array, from cell...
    if iscell(filenames_or_V), filenames_or_V = char(filenames_or_V{:}); end
    
    for i = 1:size(filenames_or_V, 1)
        currentfilename = deblank(filenames_or_V(i, :));

        if ~(exist(currentfilename, 'file'))

            fprintf('Creating Image: %3.0f ', i);
            
            if isempty(varargin) && ~exist('exampleV', 'var')
                error('Must pass in example V struct!');
            elseif ~exist('exampleV', 'var')
                exampleV = varargin{1};
            end
    
            % Create empty float image with no spm scaling (pinfo intercept is 0,
            % slope is 1)
            % --------------------------------------------------------------
            %       V.fname - the filename of the image.
            %       V.dim   - the x, y and z dimensions of the volume
            %       V.dt    - A 1x2 array.  First element is datatype (see spm_type).
            %                 The second is 1 or 0 depending on the
            %                 endian-ness. 
            %                   * 1 is big-endian, 0 is little-endian
            %       V.mat   - a 4x4 affine transformation matrix mapping from
            %                 voxel coordinates to real world coordinates.
            %       V.pinfo - plane info for each plane of the volume.
            %              V.pinfo(1,:) - scale for each plane
            %              V.pinfo(2,:) - offset for each plane
            %                 The true voxel intensities of the jth image are given
            %                 by: val*V.pinfo(1,j) + V.pinfo(2,j)
            %              V.pinfo(3,:) - offset into image (in bytes).
            %                 If the size of pinfo is 3x1, then the volume is assumed
            %                 to be contiguous and each plane has the same scalefactor
            %                 and offset.

            Vout = struct('fname', currentfilename, 'mat', exampleV.mat, 'dim', exampleV.dim, 'pinfo', [1 0 0]');
            %if isfield(exampleV, 'dt'), Vout.dt = exampleV.dt; end      % SPM5 only; data type
            if isfield(exampleV, 'n'), Vout.n = exampleV.n;  end        % SPM5 only
            
            % Make sure datatype is float / float32
            switch(spm('Ver'))
                case 'SPM2'
                    Vout.dim(4) = spm_type('float'); 
                    
                case 'SPM5'
                    Vout.dt = [spm_type('float32') 0];
                    if isfield(exampleV, 'dt'), Vout.dt(2) = exampleV.dt(2); end      % use endian-ness of example.
                    
                case 'SPM8'
                    Vout.dt = [spm_type('float32') 0];
                    if isfield(exampleV, 'dt'), Vout.dt(2) = exampleV.dt(2); end      % use endian-ness of example.
                    
                otherwise
                    error('Unknown SPM version.');
            end
            
            try
                V(i) = spm_create_vol(Vout);
                create_empty(V(i));
                
            catch
                try
                    % You cannot create only some volumes in a list, if
                    % others exist...work-around...
                    Vtmp = spm_create_vol(Vout);
                    Vtmp = spm_vol(Vtmp);
                    create_empty(Vtmp);
                    V(i) = spm_vol(Vtmp.fname);
                catch

                    disp('Cannot create vol.  This could occur for several reasons:')
                    disp('1) you may not have write permission.')
                    disp('2) Other reasons...?');
                    error('Quitting.');
                end
            end

        else
            % it exists; get it
            V(i) = spm_vol(currentfilename);
        end
    end

    
    
else
    % fprintf('Checking spm_vol structures ');
    % Make sure spm_vol structure (SPM5 only)
    % alter data type if needed
    switch(spm('Ver'))
        case 'SPM2'
            V = filenames_or_V; % Skip because hard to check.
        case {'SPM5', 'SPM8'}
            V = spm_vol(filenames_or_V);
%         otherwise 
%             error('Fix or verify image writing for SPMxx. Not implemented yet.');
        otherwise
            error('Unknown SPM version.');
    end
end

str = sprintf('Writing slice data: %04d', 0);
fprintf('%s', str)

% V is a structure array, one structure per filename
% For each file to write to...write a plane

for i = 1:length(V)
    if mod(i, 10) == 0, fprintf('\b\b\b\b%04d', i); end
    
    spm_write_plane(V(i), dat(:, :, i), wh_slice);
end

erase_string(str)
%fprintf('\n');

end



function create_empty(V)
% Write empty data
switch(spm('Ver'))
    case 'SPM2'
        spm_write_vol(V, zeros(V.dim(1:3)));
    case {'SPM5', 'SPM8'}
        % assign directly
        for i = 1:length(V), V(i).private.dat(:, :, :) = 0; end
    otherwise
        error('Unknown SPM version.')
end
end

