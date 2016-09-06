function ivecobj  = region2imagevec(cl, varargin)
% Convert a region object to an image_vector object, replacing the voxels
% and reconstructing as much info as possible. Optional: Resample to the
% space of another image_vector object.
%
% The .dat field of the new "ivecobj" is made from the cl.all_data field.
% if this is empty, uses cl.val field, then cl.Z as a backup.
% Mask information is available in ivecobj.volInfo.
%
% :Usage:
% ::
%
%    ivecobj = region2imagevec(cl, [image_vector object to resample space to])
%

ivecobj = image_vector;
ivecobj.volInfo.mat = cl(1).M;
ivecobj.volInfo.dim = cl(1).dim;

% Data in image_vector objects is stored in standard matlab vectorization order
% So cannot assume that values in all_data and XYZ match in order. Must rebuild.

% Convert to 3-d mask
% -------------------------------------------------------
[~, mask] = clusters2mask2011(cl, cl(1).dim); % 2nd output: Z-field values stored in mask elements

% Vectorize
% -------------------------------------------------------

maskvec = mask(:);

valid_vox = maskvec ~= 0 & ~isnan(maskvec);
wh_valid_vox = find(valid_vox);
n = sum(valid_vox);
  
ivecobj.dat = maskvec(valid_vox);
ivecobj.removed_voxels = ~valid_vox;

% Add all_data instead of Z if we have it, or .val or .Z
% -------------------------------------------------------
% *this bit still needs to be tested for bugs*
% alldat = cat(2, cl.all_data)';
% if ~isempty(alldat)
%     XYZ = cat(2, cl.XYZ);
%     ind = sub2ind(cl(1).dim', XYZ(1, :), XYZ(2, :), XYZ(3, :));
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
if ~strcmp(cl(1).source_images,'') && ~iscell(cl(1).source_images)
    ivecobj.volInfo.fname = ['Reconstructed from region object, source: ' cl(1).source_images(1, :)];
else
    ivecobj.volInfo.fname = 'Reconstructed from region object';
end

% resample, if asked for
% -------------------------------------------------------
if length(varargin) > 0
    
    sampleto = varargin{1};
    
    if ~isa(sampleto, 'image_vector')
        error('2nd argument must be an image_vector object (including fmri_data / statistic_image) to sample to.');
    end
    
    ivecobj = resample_space(ivecobj, sampleto);
    
end

end
