% annotate_continuous_neuroimage_maps Annotate an fmri_data object against continuous reference maps.
%
% :Usage:
% ::
%
%     % This file is currently a script that operates on a variable
%     % named `obj` already in the workspace; it is not yet wrapped as
%     % a callable method. Run it from a context where `obj` is defined:
%     obj = load_image_set('emotionreg');
%     annotate_continuous_neuroimage_maps;
%
% Annotates each image in obj by computing spatial correlations
% with a curated set of continuous-valued reference maps and a few
% Neurosynth feature batteries. The reference batteries currently
% covered are:
%
%   - Margulies principal gradient of functional connectivity
%     (transmodal vs. unimodal axis).
%   - Allen Brain transcriptomic gradients.
%   - Hansen Neuromaps PET neurochemical tracer maps.
%   - Neurosynth topics (forward and reverse inference) and term maps
%     (reverse inference).
%   - Yeo / Buckner resting-state network similarity, plotted against
%     the Yeo / Schaefer 17-network atlas.
%
% :Inputs:
%
%   **obj:**
%        fmri_data object in the caller workspace whose images are to
%        be annotated. Must be present before this script is run.
%
% :Outputs:
%
%   **t:**
%        MATLAB table accumulated in the caller workspace
%        containing one row per image and one column per annotation
%        battery (transmodal gradient values, transcriptomic PCs,
%        neurochemical tracer correlations, top Neurosynth
%        topics / terms, etc.).
%
%   Several figures are also created (Yeo / Buckner network plot,
%   wedge plot by Yeo17 networks).
%
% :Subfunctions:
%
%   - prep_and_extract_values_continuous_reference_maps resamples a
%     reference map (or stack of reference maps) to the test object's
%     space, z-scores image-wise, and returns per-image correlations
%     with the test object.
%   - plot_colored_bars renders a horizontal bar plot with one
%     colour per network from a vector of values and labels.
%
% :Notes:
%
%   - There are different procedures for smooth maps (which can be
%     resampled to the test object) and atlases (which can have very
%     small values and are high-resolution, so it is better to resample
%     the data images to them).
%   - This file does not currently begin with a function declaration
%     and therefore behaves as a MATLAB script rather than a callable
%     method on the fmri_data class. The doc block above describes the
%     intended usage; converting to a proper method is a separate task.
%
% :See also:
%   - neurosynth_feature_labels, neurosynth_lexical_plot
%   - image_similarity_plot, image_similarity_plot_bucknermaps
%   - wedge_plot_by_atlas

%% Transmodal vs. unimodal:  Principal gradient of functional connectivity

disp('Transmodal vs. unimodal:  Principal gradient of functional connectivity')

pringrad1 = load_image_set('principalgradient', 'noverbose');
pringrad1.dat = -pringrad1.dat;      % flip values so transmodal is positive here

transmodal_pringrad1 = prep_and_extract_values_continuous_reference_maps(pringrad1, obj);

% obj_resamp = resample_space(obj, pringrad1);
% transmodal_pringrad1 = canlab_pattern_similarity(obj_resamp.dat, pringrad1.dat, 'correlation');

t = table(transmodal_pringrad1);  %, 'Rownames', obj.image_names)

%% Allen brain project transcriptomic gradients

[transcriptomic_grads, transcriptomic_names] = load_image_set('transcriptomic_gradients', 'noverbose');

gene_pc_values = prep_and_extract_values_continuous_reference_maps(transcriptomic_grads, obj);

t = addvars(t, gene_pc_values);

%% Hansen Neuromaps PET neurochemical tracer maps

disp('Neuromaps PET neurochemical tracer maps')

[neurochem, neurochemnames] = load_image_set('hansen22', 'noverbose');

neurochem_vals = prep_and_extract_values_continuous_reference_maps(neurochem, obj);

% name these individually and put in different columns

t = addvars(t, neurochem_vals);



%% Neurosynth topic and term maps

disp('Neurosynth topic and term maps')

ns = load(which('neurosynth_data_obj.mat'));

%        'neurosynth_topics': 
%                   54 topic maps from Yarkoni & Poldrack 2014 topic modeling analysis
%                   Selected from 100 topics for psychological relevance
%                   and given ChatGPT-based summary topic labels by Ke et al. 2024, Nat Neurosci

disp('Neurosynth topics - forward inference maps')
[image_by_feature_correlations, top_feature_tables, top_ns_maps, bottom_ns_maps, all_ns_r_table] = neurosynth_feature_labels(obj, 'topics_fi', 'images_are_replicates', false, 'noverbose');

% collect these
terms = {};
for i = 1:size(obj.dat, 2)
    terms(i, :) = top_feature_tables{i}.Term_or_Topic_highest';
end
Neurosynth_FI_topics = terms; 
t = addvars(t, Neurosynth_FI_topics);

disp('Neurosynth topics - reverse inference maps')
[image_by_feature_correlations, top_feature_tables, top_ns_maps, bottom_ns_maps, all_ns_r_table] = neurosynth_feature_labels(obj, 'topics_ri', 'images_are_replicates', false, 'noverbose');

% collect these
terms = {};
for i = 1:size(obj.dat, 2)
    terms(i, :) = top_feature_tables{i}.Term_or_Topic_highest';
end
Neurosynth_RI_topics = terms; 
t = addvars(t, Neurosynth_RI_topics);

disp('Neurosynth terms- reverse inference maps')
[image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(obj, 'images_are_replicates', false, 'noverbose');

% collect these
terms = {};
for i = 1:size(obj.dat, 2)
    terms(i, :) = top_feature_tables{i}.Term_or_Topic_highest';
end
Neurosynth_RI_terms = terms; 
t = addvars(t, Neurosynth_RI_terms);

%% Yeo/Buckner lab 

disp('Correlations with Yeo/Bucker resting-state networks')

create_figure('Yeo-Buckner networks', 1, 2);
stats = image_similarity_plot(obj, 'bucknerlab_wholebrain', 'bicolor');

subplot(1, 2, 2)
plot_colored_bars(stats.r, stats.networknames);
set(gca, 'Position', [0.6162    0.300    0.2888    0.5150])

drawnow; snapnow;

% Yeo/Schafer 
[hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(test_mapz, 'atlases', {'yeo17networks'}, 'colors', seaborn_colors(34));

drawnow; snapnow;

% end  % main function


%%

function sim_values = prep_and_extract_values_continuous_reference_maps(reference_map, obj)
% reference_map can be one or more continuous-valued maps with annotations - e.g., principal gradients, topics, neurotransmitters
% test_map is a single binarized map to annotate

% This removes empty/0 voxels so will not work for binary reference maps.
% Different procedures are appropriate for binary reference images
numUniqueValues = size(unique(reference_map.dat(:)), 1);
if numUniqueValues <= 2
    error('Reference maps do not have enough unique values. Note: This function is not appropriate for binary maps.')
end

% Z-score reference image values to put on standard scale
% If there are multiple reference images, do this image-wise 
% This removes empty/0 voxels so will not work for binary reference maps.
% this handles excluding missing 0 values column-wise, so no need to
% remove_empty first.
reference_mapz = rescale(reference_map, 'zscoreimages');

% make sure voxel lists line up and we haven't removed empty voxels
reference_mapz = replace_empty(reference_mapz);
obj = replace_empty(obj);

% map the gradient to the test map
% (rather than the reverse, as gradient map is typically less precise
% anyway and interpolation errors matter less)

reference_mapz = resample_space(reference_mapz, obj);

sim_values = canlab_pattern_similarity(obj.dat, reference_mapz.dat, 'correlation');


end % function




function plot_colored_bars(values, names, varargin)

colors = seaborn_colors(size(values, 1));

if ~isempty(varargin)
    colors = varargin{1};
end

hBar = barh(values, 'FaceColor', 'flat');

% Apply a different color to each bar using the colors cell array
for i = 1:length(values)
    hBar.CData(i, :) = colors{i};
end

set(gca, 'YTick', 1:length(names), 'YTickLabel', names, 'FontSize', 18);
axis tight;
xlabel('Correlation')

end % plot bars function


