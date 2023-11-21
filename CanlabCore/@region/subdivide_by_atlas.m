function [r, percent_coverage] = subdivide_by_atlas(r, varargin)
% Subdivide a set of blobs stored in a region object by an anatomical atlas.
%  It takes each contiguous blob in the original region object and subdivides 
% it into a new region object with separate regions for each atlas parcel covered 
% within each blob. Thus, if your original region object has 3 blobs, and the blobs 
% span 1 5 and 2 atlas parcels, respectively, then your new region object will have 
% 8 regions. 
% - The atlas parcel labels will be included in the output r.shorttitle. 
% - Two or more parcels can have the same anatomical label. 
% - Voxels in original blobs (input r) with no anatomical labels will be excluded
%
% :Usage:
% ::
%
%    [r, percent_coverage] = subdivide_by_atlas(r, [atlas name])
%
% :Inputs:
%
%   **r:**
%        a region object, defined using region(mask)
%
%   **atlas name:**
%        atlas-class object to subdivide by; see load_atlas() for options.
%
% :Output:
%
%   **r:**
%   A region object with separate clusters for each contiguous blob,
%   subdivided by regions labeled in atlas.
%
%   **percent_coverage:**
%   Percentage of voxels in the atlas parcel covered by the region
%
% :Example:
% ::
%
%    r = subdivide_by_atlas(r);
%    r(cat(1, r.numVox) < 20) = []; % get rid of small regions
%    cluster_orthviews(r, 'unique');
%

% Programmers notes:
% Changed 4/23/21 by Tor Wager
% To take atlas object as input, NOT file.


percent_coverage = [];

if isempty(r), return, end  % for case when calling with empty clusters 

if ~isa(r, 'region'), error('Object is not a region object'); end

% The cases above shouldn't happen, but they apparently can when calling from some
% enclosing functions

if nargin < 2 || isempty(varargin{1})
    error('Enter an atlas object as the 2nd input. See load_atlas( ) for some options.');
    % atlasname = which('atlas_labels_combined.img');
else
    atlas_obj = varargin{1};
end

%     if ~exist(atlasname, 'file')
%         if exist(which(atlasname), 'file')
%             atlasname = which(atlasname);
%         else
%             fprintf('Cannot find atlas image on path! Looking for:%s\n', atlasname);
%         end
%     end

%r = region(overlapmask); % r is input

% if is region, reconstruct...
% ivec is region object voxels (the full set) reconstructed into image_vector
[ivec, orig_indx] = region2imagevec(r);

% region2fmri_data

% region2imagevec creates illegal list of removed_voxels (doesn't match
% .volInfo.wh_inmask).  Fix...
% ivec.removed_voxels = ivec.removed_voxels(ivec.volInfo.wh_inmask);

% label_mask = fmri_data(atlasname);

% resample and mask label image
original_atlas_obj = atlas_obj;                         % save for calculating coverage
atlas_obj = resample_space(atlas_obj, ivec, 'nearest'); % this eliminates atlas regions

ulabels = unique(atlas_obj.dat(atlas_obj.dat ~= 0));

fprintf('%3.0f unique non-zero atlas labels covered by the input region object (label = 0 will be excluded).\n', length(ulabels));
fprintf('%3.0f voxels in the input region object have atlas labels of 0 and will not be labeled.\n', length(atlas_obj.dat == 0));

% Define regions based on unique voxels values in mask_image, and extract
% data stored in data_comb object, resampled to space of mask_image.
% Space is defined by mask_image:\

% for each of the original regions...
for i = 1:length(r)
    %ivec = region2imagevec(r(i));
    %ivec = remove_empty(ivec);
    
    % orig_indx is the original cluster ID
    whvox = orig_indx == i;
    nvox = sum(whvox);
    
    if nvox == 1
        rr{i} = r(i); % just copy - nothing to subdivide
        
        % add title
        %wh = atlas_obj.dat(whvox);
        label_indx = atlas_obj.dat(orig_indx == i);
        if label_indx == 0
            rr{i} = [];
        else
            rr{i}.shorttitle = atlas_obj.labels{label_indx};
        end
        
    else
        %mylabel = label_mask;
        mylabel = atlas_obj;
        mylabel.dat(~whvox, 1) = 0; % remove these - not in contiguous region
        
        if ~any(mylabel.dat)
            % no labels - exclude
            rr{i} = [];
        else
            
            % Get all atlas regions covered by the original region, return
            % regions in rr{i} with names in rr{i}.shorttitle
            
            % mylabel is atlas object with .dat containing new anatomical labels
            % ivec is image_vector object with voxels in input region blobs
            % voxel lists/spaces match
            
            rr{i} = atlas2region(mylabel, ivec, 'unique_mask_values', 'noverbose');
            
            for jj = 1:length(rr{i})
                if size(rr{i}(jj).XYZ, 2) ~= size(rr{i}(jj).val, 1)
                    % illegal
                    warning('Illegal length for data values...do not match voxels.')
                end
            end
            
        end
        
    end
    
end

% Count how many original regions do not have any labeled voxels 
for i = 1:length(rr)
    wh_omit(1, i) = isempty(rr{i});
end

fprintf('%3.0f original regions do not have any labeled voxels.\n', sum(wh_omit));

r = cat(2, rr{:});

% If requested, calculate coverage of each labeled region -- percentage of
% voxels in the atlas parcel in the region

if nargout > 1

    percent_coverage = zeros(length(r), 1);

    % relative volume of a voxel in atlas compared to region
    conversion_factor = prod(diag(abs(original_atlas_obj.volInfo.mat))) ./ prod(diag(abs(r(1).M)));

    for i = 1:length(r)

        parcel = select_atlas_subset(original_atlas_obj, {r(i).shorttitle});

        percent_coverage(i) = 100 .* r(i).numVox ./ (sum(parcel.dat > 0) .* conversion_factor);

    end

end


end % function

