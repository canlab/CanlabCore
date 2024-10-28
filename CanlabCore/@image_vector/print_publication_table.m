function [outputT_pos, outputT_neg] = print_publication_table(image_vec, varargin)
% This function is designed to generate an fMRI activation table for use in
% main texts of papers.
% 
% Specifically, the tables generated have the following characteristics
% (1) Each row of the table is one cluster, or "blob", of activations
% (2) Columns: atlas labels, heuristic names (now only work for canlab2023), 
% peak coordinates, and t values
% (3) Atlas labels and heuristic names can have multiple values, separate
% by commas
%
% This function works by concatenating the ouputs from
% `autolabel_regions_using_atlas()` and `table()`
% 
% Zizhuang Miao, Michael Sun, 09/17/2024
%
% :Usage:
% ::
%
%    [outputT_pos, outputT_neg] = print_publication_table(image_vec, varargin)
%
% :Inputs:
%
%    **image_vec**:
%         An image vector object
%
% :Optional keyword arguments:
%
%    **'atlas'**:
%        Followed by a char array with the keyword or name of an 
%        atlas-class object with labels.
%        If not provided, default to 'canlab2023'
%
%    **'doneg'**:
%        Return both tables with positive activations and tables with 
%        negative activations. Default to false.
%    
%    **'noverbose'**:
%        Suppress the output tables (only work on the modified `table()`
%        function)
%
%    **'nolegend'**:
%        Suppress table legends
%
%    **'threshold'**:
%        Followed by the threshold for the atlas. Default to 0.
%
% :Outputs:
%
%    **outputT_pos**
%        A Matlab table object with one row being a cluster/blob of
%        positive activation.
% 
%    **outputT_neg**
%        Same table as above, with negative activation.
% ..
%    Programmers' notes:
%
%    10/28/2024  Zizhuang Miao: Added warnings when the space of the image
%    vector does not match that of the atlas. It could result in more false
%    positives in the output table. It is rooted in the way
%    `autolabel_regions_using_atlas` works, which calls `resample_space`
%    and `probability_maps_to_region_index` subsequently.
%    'probability_maps_to_region_index' works in a winner-take-all way such
%    that as long as a voxel is close to a significant voxel, it can 
%    become parts of a significant region in the output table. It is not
%    preferrable and should be avoided. Now a heuristic way is to resample
%    the image vector as soon as possible (e.g., before running t tests).

narginchk(1, Inf)


if ~isa(image_vec, 'region')

    if (size(image_vec.dat, 2) > 1)
        error('Use this function with fmri_data or image_vector objects containing only a single image. Use get_wh_image() to select one.')
    end

end
% Default
atlas_obj = [];
nTables = 1;    % the number of tables to produce
                % 1 for positive only, 2 for both positive and negative
verbosestr = "";
legendstr = "";
thr = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'atlas'}, atlas_obj = load_atlas(varargin{i+1});
            
            case {'doneg'}, nTables = 2;
            
            case {'noverbose'}, verbosestr = ", 'noverbose'";

            case {'nolegend'}, legendstr = ", 'nolegend'";

            case {'threshold'}, thr = varargin{i+1};
        end
    end
end

% -------------------------------------------------------------------------
% PREP ATLAS AND OBJECT
% -------------------------------------------------------------------------

% Load default atlas if needed
if isempty(atlas_obj)
    % atlas_obj = load_atlas('canlab2023');
    atlas_obj = load_atlas('canlab2024');
end

% threshold atlas
atlas_thr = atlas_obj.threshold(thr);

% find network labels for all regions
% iterate over labels and find corresponding network labels from canlab2018
% cr. Michael Sun

disp('Generating Network Column ...');

canlab2018 = load_atlas('canlab2018');
nets = canlab2018.labels_2;
labels = atlas_thr.labels;
sorted_nets = cell(1,length(labels));
for i = 1:length(labels)
    if contains(labels{i},'Ctx')
        sorted_nets{i} = nets{strcmp(canlab2018.labels,labels{i})};
    else
        sorted_nets{i} = 'Sub-cortex';
    end
end
atlas_thr.labels_5 = sorted_nets;    % labels_5 should be a safe field name
                                     % with no values in the field before
% create a new atlas object with the network labels
atlas_net = atlas_thr.downsample_parcellation('labels_5');    

% parcellate image_vec into (positve and negative) regions
if ~isa(image_vec, 'region')
    regions = region(image_vec);
else
    regions = image_vec;
end

% check whether the space of the regions and atlas match
if any(regions.volInfo.dim ~= atlas_thr.volInfo.dim) || ...
        any(regionsvolInfo.mat ~= atlas_thr.volInfo.mat)
    warning(['The space of the image vector does not match that of the atlas. ' ...
        'This can result in false positives in output tables. ' ...
        'We recommend resampling the image to the atlas ' ...
        'space as early as possible, not at this stage.']);
end

r = cell(1, 2);    % a cell object that stores positive and negative regions
r_relabeled = cell(1, 2);
outputTables = cell(1, 2);   % a cell that stores output tables
[r{1}, r{2}] = posneg_separate(regions);

% -------------------------------------------------------------------------
% CREATE AND CONCATENATE TABLES
% -------------------------------------------------------------------------
for t = 1:nTables
    % first, use `autolabel_...` to get each blob and names of atlas regions
    [r_relabeled{t}, atlasT, ~, atlasR] ...
        = autolabel_regions_using_atlas(r{t}, atlas_thr);

    % only retain blobs that covered at least one atlas region
    atlasR = atlasR(atlasT.Atlas_regions_covered>0, :);
    atlasT = atlasT(atlasT.Atlas_regions_covered>0, :);
    
    % next, create a table with network names to be attached to the final table
    [~, networkT] = autolabel_regions_using_atlas(r_relabeled{t}, atlas_net);
    
    % then, find peak coordinates and stats values for each blob
    eval(strcat("[~, ~, xyzT] = table(r_relabeled{t}, 'atlas_obj', atlas_thr, 'noverbose'", ...
        legendstr, verbosestr, ");"))

    atlasT.heuristic_names = cell(height(atlasT), 1);

    % add peak coordinates and atlas region names
    [atlasT.x, atlasT.y, atlasT.z, atlasT.Region, atlasT.t] ...
        = deal(NaN(height(atlasT), 1));
    for i = 1:height(atlasT)
        atlasT.Region(i) = i;     % name the regions by number
        atlasT.all_regions_covered{i} = atlasR{i};
        % find corresponding clusters across two tables by volumes and #regions
        nVolume = atlasT.Region_Vol_mm(i);
        nRegions = atlasT.Atlas_regions_covered(i);
        perc = atlasT.Perc_covered_by_label(i);
        wh = find(xyzT.Volume == nVolume & ...
            xyzT.Atlas_regions_covered == nRegions & ...
            xyzT.Perc_covered_by_label == perc);
        atlasT.x(i) = xyzT.XYZ(wh, 1);
        atlasT.y(i) = xyzT.XYZ(wh, 2);
        atlasT.z(i) = xyzT.XYZ(wh, 3);
        atlasT.t(i) = xyzT.maxZ(i);
    end
    
    % add network names and heuristic names (now only work for canlab2023)
    atlasT.network = cell(height(atlasT), 1);
    for i = 1:height(atlasT)
        % network
        nVolume = atlasT.Region_Vol_mm(i);
        nRegions = atlasT.Atlas_regions_covered(i);
        wh_networks = find(networkT.Region_Vol_mm == nVolume);
        atlasT.network{i} = networkT.modal_label{wh_networks};
        
        % heuristic names (only work for 'canlab2023' now)
        % if strcmp(atlas_obj.atlas_name, 'CANLab2023_MNI152NLin2009cAsym_coarse_2mm')
            heuName = cell(1, size(atlasR{i}, 2));
            for j = 1:size(atlasR{i}, 2)
                heuName{j} = atlas_thr.labels_3{strcmp(atlas_thr.labels, atlasR{i}{j})};
            end
            heuName = unique(heuName, 'stable');
            atlasT.heuristic_names{i} = heuName;
        % end
    end
    
    % remove useless columns
    atlasT = removevars(atlasT, ...
        {'Voxels', 'Atlas_regions_covered', 'modal_label', ...
        'modal_label_descriptions', 'Perc_covered_by_label', ...
        'Ref_region_perc', 'modal_atlas_index'});
   
    % remove all underlines in region names, and put them as a long string
    atlasT.all_regions_covered = cellfun(@(x) strjoin(x, ', '), ...
        format_strings_for_legend(atlasT.all_regions_covered), ...
        'UniformOutput', false);
    atlasT.heuristic_names = cellfun(@(x) strjoin(x, ', '), ...
        format_strings_for_legend(atlasT.heuristic_names), ...
        'UniformOutput', false);

    % rename columns
    atlasT.Properties.VariableNames = ...
        {'Cluster', 'Volume (mm^3)', 'Atlas region names', ...
        'Heuristic names', 'X', 'Y', 'Z', 'Max t', 'Network'};
    
    if isa(image_vec, 'statistic_image')
        if image_vec.type ~= 'T'
            atlasT = renamevars(atlasT, {'Max t'}, {['Max ', lower(image_vec.type)]});
        end
    end
    outputTables{t} = atlasT;
end

outputT_pos = outputTables{1};
outputT_neg = outputTables{2};

end