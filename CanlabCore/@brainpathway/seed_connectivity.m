% Estimate voxel-wise connectivity with one or more 'seed' regions
%
% :Usage:
% ::
%
% fmri_dat_connectivity_maps = seed_connectivity(obj, varargin)
%
% Given a brainpathway object, return voxel-wise maps of connectivity with
% one or more 'seed' regions defined in the region_atlas atlas object.
% - Uses Pearson's correlation metric on voxel data (.voxel_dat field)
% - Returns output in an fmri_data object
% - By default, returns connectivity maps for ALL regions in the atlas
%
% Dev notes:
% - this needs more development to work for nodes if nodes/regions do not match;see find_node_indices
% - uses Pearson's correlations now; needs more dev to use metric stored in object.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2019 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        A brainpathway object
%
% :Optional Inputs:
%   **MULTIPLE:**
%        Inputs to atlas.select_atlas_subset
%        Generally, this will be a cell array with one or more names of
%        regions or nodes (must match those in obj.region_atlas.labels and/or
%        obj.node_labels) ... or a vector of integers indicating which
%        regions/nodes.
%
%   **'regions':**
%       Calculate connectivity among region averages only (default)
%
%   **'nodes':**
%       Calculate connectivity among node responses only
%
%   **param2:**
%        future: could change correlation metrics or other options
%
% :Outputs:
%
%   **fmri_dat_connectivity_maps:**
%        An fmri_data object with one image per seed region selected
%        .dat field contains correlation values
%
%   **rr:**
%       Raw correlation values, not expanded (e.g., regions x images).
%
% :Examples:
% ::
%
% Assume you have an fmri_data object with time series data in
% img_resampled, and that this is in the same space as the brainpathway
% region_atlas object (see resample_space method).
%
% % 1. Define a brainpathway object 'regions' consisting of 32 networks
%    b = brainpathway(load_atlas('yeo17networks'));
%
% % 2. Attach the data to extract
%    b.voxel_dat = img_resampled.dat;
%
% % 3. Calculate correlation maps with all network averages containing
% % 'Default' in the label
%   fmri_dat_connectivity_maps = seed_connectivity(b, {'Default'});
%
% % 4. Visualize the results
%   orthviews(fmri_dat_connectivity_maps);
%
% % 5. Now replace the atlas with the CANlab 500-region atlas, and estimate
% % connectivity maps for left and right Nuc Accumbens (NAC):
%   b.region_atlas = load_atlas('canlab2018_2mm');
%   fmri_dat_connectivity_maps = seed_connectivity(b, {'NAC'});
%   orthviews(fmri_dat_connectivity_maps);
%
% :References:
%   None needed.
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

function [fmri_dat_connectivity_maps, rr] = seed_connectivity(obj, varargin)

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

regions_or_nodes = 'regions';

% This function takes and N x p matrix a and an N x v matrix b and returns
% a p x v matrix of correlations across the pairs.
corr_matrix = @(a, b) ((a-mean(a))' * (b-mean(b)) ./ (size(a, 1) - 1)) ./ (std(b)' * std(a))'; % Correlation of a with each column of b

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'region', 'regions'}, regions_or_nodes = 'regions'; varargin{i+1} = [];
            case {'node', 'nodes'}, regions_or_nodes = 'nodes'; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% Pass all inputs to select_atlas_subset
% NOTE: this needs more development to work for nodes if nodes/regions do
% not match ****** see find_node_indices
[~, to_extract] = select_atlas_subset(obj.region_atlas, varargin{:});

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

switch regions_or_nodes
    
        case 'voxels'
        % ------------------------------------------
        
        
        if obj.verbose, fprintf('Calculating correlations bewtween voxels and seed regions.\n'); end
        
        
        a = double(obj.region_dat(:, to_extract));
        b = double(obj.voxel_dat');
        rr = corr_matrix(a, b)';    % Voxels x regions
        
        fmri_dat_connectivity_maps = fmri_data(image_vector('volInfo', obj.region_atlas.volInfo, 'dat', rr));
        
        fmri_dat_connectivity_maps.image_names = char(obj.region_atlas.labels{to_extract});
        fmri_dat_connectivity_maps.source_notes = 'Voxel-wise correlation maps for seed region averages created with brainpathway.seed_connectivity';
        
        
    case 'regions'
        % ------------------------------------------
        
        
        if obj.verbose, fprintf('Calculating correlations between region averages and seed regions.\n'); end
        
        a = double(obj.region_dat(:, to_extract));
        b = double(obj.region_dat);
        rr = corr_matrix(a, b)';    % Voxels x regions
        
        % Need to expand to correlation value for each voxel
        rr_expanded = expand_values_region2voxel(obj, rr);
        
        fmri_dat_connectivity_maps = fmri_data(image_vector('volInfo', obj.region_atlas.volInfo, 'dat', rr_expanded));
        
        fmri_dat_connectivity_maps.image_names = char(obj.region_atlas.labels{to_extract});
        fmri_dat_connectivity_maps.source_notes = 'Region-wise correlation maps for seed region averages created with brainpathway.seed_connectivity';
        
    case 'nodes'
        % ------------------------------------------
        
        if isempty(obj.node_dat)
            fmri_dat_connectivity_maps = [];
            if obj.verbose, fprintf('No node data found. Skipping seed correlations with nodes.\n'); end
            
            return
        else
            if obj.verbose, fprintf('Calculating correlations between node responses and seed nodes.\n'); end
        end
        
        to_extract = find_node_indices(obj, varargin{:});
        
        a = double(obj.node_dat(:, to_extract));
        
        rr = corr_matrix(a, b)';    % Voxels x regions
        
        fmri_dat_connectivity_maps = fmri_data(image_vector('volInfo', obj.region_atlas.volInfo, 'dat', rr));
        
        fmri_dat_connectivity_maps.image_names = char(obj.node_labels{to_extract});
        fmri_dat_connectivity_maps.source_notes = 'Node-wise correlation maps for node responses created with brainpathway.seed_connectivity';
        
end % switch

end % main function


function to_extract = find_node_indices(obj, varargin)
% Identify nodes to extract from brainpathway object given an integer
% vector of which nodes to select, cell array of strings for node names,
% or combination of the two. Returns a logical vector of nodes in set
% 
% obj : a brainpathway object
% 
% Optional inputs:
% a cell array of strings containing node names (matched to obj.node_labels)
% an integer vector of which nodes
%
% 'flatten' (not used here - would collapse atlas regions if used in select_atlas_regions)

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

strings_to_find = [];
integers_to_find = [];
doflatten = false;

% optional inputs with default values
for i = 1:length(varargin)
    
    if iscell(varargin{i})
        strings_to_find = varargin{i};
        
    elseif isnumeric(varargin{i})
        integers_to_find = varargin{i};
        
    elseif ischar(varargin{i})
        switch varargin{i}
            
            case 'flatten', doflatten = true;
                
                %             case 'xxx', xxx = varargin{i+1}; varargin{i+1} = [];
                
            otherwise
                warning(['Unknown input string option:' varargin{i} '. Assuming it might be an atlas label. Place atlas labels in a cell array']);
                strings_to_find{1} = varargin{i};
        end
    end
end


% -------------------------------------------------------------------------
% INIT
% -------------------------------------------------------------------------

k = size(obj.node_dat, 2);
to_extract = false(1, k);

% -------------------------------------------------------------------------
% FIND BY STRING
% -------------------------------------------------------------------------

for i = 1:length(strings_to_find)
    
    % Find which names match
    wh = ~cellfun(@isempty, strfind(obj.node_labels, strings_to_find{i}));
    
    to_extract = to_extract | wh;
    
end

% -------------------------------------------------------------------------
% FIND BY NUMBERS
% -------------------------------------------------------------------------

to_extract(integers_to_find) = true;

if ~any(to_extract)
    error('No nodes identified to extract.');
end

end % subfunction




function output_dat = expand_values_region2voxel(obj, input_dat)
% Map region-wise values to voxel values
% output_dat = expand_values_region2voxel(obj, input_dat)
%
% obj: A brainpathway object
% input_dat: regions x images values to expand
% output_dat: voxels x images
%
% Tor Wager, 12/2019

% Integer codes for which voxels belong to each region
voxel_indx = obj.region_atlas.dat;

[n_regions, n_images] = size(input_dat);

output_dat = zeros([size(voxel_indx, 1), n_images]);


for i = 1:n_regions
    
    wh = voxel_indx == i;
    
    for j = 1:n_images
        
        myvalue = input_dat(i, j);
        
        output_dat(wh, j) = myvalue;
        
    end
    
end

end % function

