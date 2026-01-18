

% there are different procedures for smooth maps, which can be resampled to the
% object, and atlases, which can have very small values and are high-res,
% so better to sample the data images to them.

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


