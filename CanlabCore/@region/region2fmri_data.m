function objout = region2fmri_data(r, reference_obj)
% Transform region object r into fmri_data object objout, given reference
% fmri_data object reference_obj in the same space as region object
%
% objout = region2fmri_data(r, reference_obj)
%
% Examples:
% see parcellate_pain_predictive_regions.m
% lindquist2015_pos_region_obj = region2fmri_data(r, obj);
% back to regions to test:
% rtest = region(lindquist2015_pos_region_obj, 'unique_mask_values');
% orthviews(rtest, 'unique')


%objout = fmri_data();

% Complete solution: rebuild without reference mask
% Here: need reference mask
%          fname: 'REMOVED: CHANGED SPACE'
%            dim: [91 109 91]
%             dt: [2 0]
%          pinfo: [3×1 double]
%            mat: [4×4 double]
%              n: [1 1]
%        descrip: 'Space of /Users/tor/Documents/Code_External/spm12/toolbox/FieldMap/bra?'
%        private: [1×1 nifti]
%           nvox: 902629
%     image_indx: [902629×1 logical]
%      wh_inmask: [352328×1 double]
%       n_inmask: 352328
%        xyzlist: [352328×3 double]
%        cluster: [352328×1 double]

% Added data mapping. .dat transferred first, otherwise, .Z, followed by
% the region number
 
if nargin < 2
    % no reference object
    % use new 2025 simple build from scratch

    objout = region2fmri_data2025(r);

    return % skip everything below

end

isbad = false;

if reference_obj.volInfo.dim ~= r(1).dim, isbad = true; end
if any(reference_obj.volInfo.mat(:) ~= r(1).M(:)), isbad = true; end

if isbad, error('reference_obj must have same mat file and dim as region_object'); end


if ~isa(reference_obj, 'fmri_data')
    % First attempt to cast the reference object into a valid fmri_data object.
    reference_obj = fmri_data(reference_obj);
end

% rebuild from reference obj
reference_obj = replace_empty(reference_obj);
reference_obj.dat = zeros(size(reference_obj.dat, 1), 1);

wh_not_in_mask = ~reference_obj.volInfo.image_indx;

for i = 1:length(r)
    
    ivec = region2imagevec2tmp(r(i));
    
    % get voxel values in whole-image space
    wh = ivec.volInfo.image_indx;
    %  wh(find(wh)) = ivec.dat;  % all 1, not needed
    
    % get in region_obj space
    wh(wh_not_in_mask) = [];
    
    if ~isempty(r(i).dat)
        reference_obj.dat(find(wh), 1) = r(i).dat'; % code with dat
    elseif ~isempty(r(i).Z)
        reference_obj.dat(find(wh), 1) = r(i).Z'; % code with Z
    else
        reference_obj.dat(find(wh), 1) = i; % code with region number
    end
end

objout = reference_obj;

end % main function




function out_obj = region2fmri_data2025(r)

dim = r(1).dim;
nvox = prod(dim);

vol = zeros(dim);

xyz = cat(2, r.XYZ);

% remove redundant xyz
[xyz, indx] = unique(xyz', 'rows', 'stable');
xyz = xyz';

% Z must be row vector or this will not work
% note: if there are redundant voxels xyz in regions(k), this function will pick the .Z field value
% for the first one
vals = cat(2, r.Z);
vals = vals(indx);

ind = sub2ind(dim, xyz(1, :)', xyz(2, :)', xyz(3, :)');

mat = r(1).M;

vol(ind) = vals;
vol_vec = vol(:);           % valued
mask_indx = vol_vec ~= 0;   % logical

% Define volinfo
% -------------------------------------

volInfo = struct();

volInfo.fname = 'Created from region object';
volInfo.dim = dim;
volInfo.dt = [2 0];
volInfo.pinfo = [1 0 352]';
volInfo.mat = mat;
volInfo.n = [1 1];
volInfo.descrip = 'Created from region object';
volInfo.private = [];
volInfo.nvox = nvox;
volInfo.image_indx = mask_indx; %true(nvox, 1);
volInfo.wh_inmask = find(mask_indx);
volInfo.n_inmask = length(volInfo.wh_inmask);
volInfo.xyzlist = xyz';
volInfo.cluster = spm_clusters(xyz);

% define mask
% -------------------------------------

mask = fmri_mask_image;
mask.dat = single(mask_indx(mask_indx) ~= 0);
% mask.volInfo = volInfo; % May need to define more, not sure

% define object
% -------------------------------------
out_obj = fmri_data;

out_obj.mask = mask;

out_obj.mask_descrip = 'Copied from region() object';
out_obj.source_notes = 'Copied from region() object';
out_obj.dat = vals'; vol_vec(mask_indx);
out_obj.volInfo = volInfo;
out_obj.removed_voxels = 0;
out_obj.removed_images = 0;
out_obj.image_names = '';
out_obj.fullpath = '';
out_obj.files_exist = 0;
out_obj.history = {'created from region object'};

end % region2fmridata2025 function
