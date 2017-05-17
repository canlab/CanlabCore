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
 
isbad = false;

if reference_obj.volInfo.dim ~= r(1).dim, isbad = true; end
if any(reference_obj.volInfo.mat(:) ~= r(1).M(:)), isbad = true; end

if isbad, error('reference_obj must have same mat file and dim as region_object'); end

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
    
    reference_obj.dat(find(wh), 1) = i; % code with region number
    
end

objout = reference_obj;

end % main function

