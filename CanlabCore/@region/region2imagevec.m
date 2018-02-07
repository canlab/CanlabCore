function [ivecobj, orig_idx_vec]  = region2imagevec(r, varargin)
% Convert a region object to an image_vector object, replacing the voxels
% and reconstructing as much info as possible. Optional: Resample to the
% space of another image_vector object specified by the first additional input.
%
% The .dat field of the new "ivecobj" is made from the r.all_data field.
% if this is empty, uses r.val field, then r.Z as a backup.
% Mask information is available in ivecobj.volInfo.
%
% :Usage:
% ::
%
%    [ivecobj, orig_idx] = region2imagevec(r, [image_vector object to resample space to])
%

% Alt code - would need to expand for all regions
% ivecobj = image_vector;
% ivecobj.volInfo.dim = r(1).dim;
% ivecobj.volInfo.mat = r(1).M;
% ivecobj.volInfo.n = [1 1];
% ivecobj.volInfo.nvox = prod(r(1).dim);
% ivecobj.volInfo.image_indx = false(ivecobj.volInfo.nvox, 1);
% ivecobj.volInfo.fname = [];
%
% % get locations in image
% 
% voldata = zeros(r(1).dim);
% voxels = r(1).XYZ';  % columns
% 
% for i = 1:size(voxels, 1)
%     voldata(voxels(i,1),voxels(i,2),voxels(i,3)) = r(1).val(i);
% end
% 
% ivecobj.volInfo.image_indx = voldata(:) ~= 0;  % logical
% ivecobj.volInfo.wh_inmask = find(ivecobj.volInfo.image_indx);
% ivecobj.volInfo.n_inmask = length(ivecobj.volInfo.wh_inmask);
% ivecobj.volInfo.xyzlist = voxels;
% ivecobj.volInfo.cluster = ones(ivecobj.volInfo.n_inmask, 1);
% ivecobj.volInfo.descrip = sprintf('From region object %s', r(1).title);

% Initialize minimum viable image vector object
% -------------------------------------------------------
ivecobj = image_vector;
ivecobj.volInfo.mat = r(1).M;
ivecobj.volInfo.dim = r(1).dim;
ivecobj.volInfo.n = [1 1];
ivecobj.volInfo.nvox = prod(r(1).dim);
ivecobj.volInfo.image_indx = false(ivecobj.volInfo.nvox, 1);
ivecobj.volInfo.fname = [];

% Data in image_vector objects is stored in standard matlab vectorization order
% So cannot assume that values in all_data and XYZ match in order. Must rebuild.

% Convert to 3-d mask
% -------------------------------------------------------
[orig_idx, mask] = clusters2mask2011(r, r(1).dim); % 2nd output: Z-field values stored in mask elements

% Vectorize
% -------------------------------------------------------

maskvec = mask(:);
orig_idx_vec = orig_idx(:);

valid_vox = maskvec ~= 0 & ~isnan(maskvec);
wh_valid_vox = find(valid_vox);
n = sum(valid_vox);
  
ivecobj.dat = maskvec(valid_vox);
ivecobj.removed_voxels = 0;             % ~valid_vox; tor edit jan 2018. valid vox should either be in volInfo or removed_, not both
orig_idx_vec = orig_idx_vec(valid_vox);

% Add all_data instead of Z if we have it, or .val or .Z
% -------------------------------------------------------
% *this bit still needs to be tested for bugs*
% alldat = cat(2, r.all_data)';
% if ~isempty(alldat)
%     XYZ = cat(2, r.XYZ);
%     ind = sub2ind(r(1).dim', XYZ(1, :), XYZ(2, :), XYZ(3, :));
% 
%     ivecobj.dat(valid_vox, :) = alldat';
% end

ivecobj.volInfo.image_indx = valid_vox;
ivecobj.volInfo.n_inmask = n;

ivecobj.volInfo.wh_inmask = wh_valid_vox;

ivecobj.volInfo.nvox = prod(ivecobj.volInfo.dim);

[i, j, k] = ind2sub(ivecobj.volInfo.dim, wh_valid_vox);
ivecobj.volInfo.xyzlist = [i j k];

% % re-build contiguity cluster indices

if ivecobj.volInfo.n_inmask < 50000
    ivecobj.volInfo.cluster = spm_clusters(ivecobj.volInfo.xyzlist')';
else
    ivecobj.volInfo.cluster = ones(ivecobj.volInfo.n_inmask, 1);
end

% for reslicing compatibility

ivecobj.volInfo.dt = [16 0]; 

% meta-data 
if ~strcmp(r(1).source_images,'') && ~iscell(r(1).source_images)
    ivecobj.volInfo.fname = ['Reconstructed from region object, source: ' r(1).source_images(1, :)];
else
    ivecobj.volInfo.fname = 'Reconstructed from region object';
end

% resample, if asked for
% -------------------------------------------------------
if ~isempty(varargin)
    
    sampleto = varargin{1};
    
    if ~isa(sampleto, 'image_vector')
        error('2nd argument must be an image_vector object (including fmri_data / statistic_image) to sample to.');
    end
    
    ivecobj = resample_space(ivecobj, sampleto);
    iveorigidx = ivecobj;
    % same for orig idx vector
    iveorigidx.dat = orig_idx_vec;
    iveorigidx = resample_space(iveorigidx, sampleto);
    orig_idx_vec = iveorigidx.dat;
end

end
