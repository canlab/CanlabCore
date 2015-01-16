function r = subdivide_by_atlas(r, varargin)
% r = subdivide_by_atlas(r, [atlas name])
%
% Inputs:
%   r = a region object, defined using region(mask)
%   atlas name: Optional mask image with integer codes defining in-mask
%   regions.  Default is 'atlas_labels_combined.img'
%
% Output: A region object with separate clusters for each contiguous blob,
% subdivided by regions labeled in atlas.
%
% Example:
% r = subdivide_by_atlas(r);
% r(cat(1, r.numVox) < 20) = []; % get rid of small regions
% cluster_orthviews(r, 'unique');

if nargin < 2 || isempty(varargin{1})
atlasname = which('atlas_labels_combined.img');
else
    atlasname = varargin{1};
end

if ~exist(atlasname, 'file')
    if exist(which(atlasname), 'file')
        atlasname = which(atlasname);
    else
        fprintf('Cannot find atlas image on path! Looking for:%s\n', atlasname);
    end
end


%r = region(overlapmask); % r is input

% if is region, reconstruct...
[ivec, orig_indx] = region2imagevec(r);

label_mask = fmri_data(atlasname);

% resample and mask label image
label_mask = resample_space(label_mask, ivec, 'nearest');

ulabels = unique(label_mask.dat);

fprintf('%3.0f unique atlas labels in mask (label = 0 will be excluded).\n', length(ulabels));

% Define regions based on unique voxels values in mask_image, and extract
% data stored in data_comb object, resampled to space of mask_image.
% Space is defined by mask_image:
for i = 1:length(r)
    %ivec = region2imagevec(r(i));
    %ivec = remove_empty(ivec);
    
    whvox = orig_indx == i;
    nvox = sum(whvox);
    
    if nvox == 1
        rr{i} = r(i); % just copy - nothing to subdivide
    else
        mylabel = label_mask;
        mylabel.dat(orig_indx ~= i, 1) = 0; % remove these - not in contiguous region
        
        if ~any(mylabel.dat)
            % no labels - exclude
            rr{i} = [];
        else
            rr{i} = region(mylabel, ivec, 'unique_mask_values');
        end
        
    end
    
end

for i = 1:length(rr)
wh_omit(1, i) = isempty(rr{i});
end

fprintf('%3.0f original regions do not have any labeled voxels.\n', sum(wh_omit));

r = cat(2, rr{:});


end % function

