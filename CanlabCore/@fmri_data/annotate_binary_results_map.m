function RESULTS = annotate_binary_results_map(test_map)
% Annotate a binary results map with several types of gradients and networks
%
% - Transmodal vs. unimodal:  Principal gradient of functional connectivity
% - Allen brain project transcriptomic gradients
% - Hansen Neuromaps PET neurochemical tracer maps
% - Neurosynth topics and terms
% - Yeo/Buckner resting-state fMRI maps
%
% :Usage:
% ::
%
%     annotate_binary_results_map(fmri_data_object)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2025 Tor Wager
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
%   **fmri_data_object:**
%        fmri_data object with a binary map
%
% :Optional Inputs:
%   **none yet:**
%
% :Outputs:
%
%   **none yet**
%
% :Examples:
% ::
%
% ns = load(which('neurosynth_data_obj.mat'));
% test_map = get_wh_image(ns.topic_obj_reverseinference, 1); % somatosensory topic
% annotate_binary_results_map(test_map)
%
% :See also:
%   image_similarity_plot, wedge_plot_by_atlas

% uses:
% kernelDensitySmoothedHistogram, plotMatrixKernelDensity


%% Prep map for correlation analyses and plots

% Note: we must z-score the test map here to prevent downstream functions
% from treating zeros in the test_map as missing values.  
% - This is important when calculating correlations with topic maps, as excluding missing
% values from the test map will result in correlations calculated only on
% the variation within non-zero values, which may not be meaningful.
% - Using rescale/zscoreimages will not work for this purpose, as it
% zscores within non-zero values, assuming 0 is a not-a-number value (SPM
% convention)

test_mapz = test_map;
test_mapz.dat = zscore(test_mapz.dat);
% test_mapz = rescale(test_map, 'zscoreimages');



%% Transmodal vs. unimodal:  Principal gradient of functional connectivity

disp('Transmodal vs. unimodal:  Principal gradient of functional connectivity')

pringrad1 = load_image_set('principalgradient', 'noverbose');
pringrad1.dat = -pringrad1.dat;      % flip values so transmodal is positive here

% surface(pringrad1, 'foursurfaces')
% surface(pringrad1, 'hcp inflated left');

% reference_map (pringrad1) can be one or more continuous-valued maps with annotations - e.g., principal gradients, topics, neurotransmitters
% test_map is a single binarized map to annotate
[vals, vals_negative, minPlotVal, maxPlotVal] = prep_and_extract_values_continuous_reference_maps(pringrad1, test_map);

[density_vals, xValues] = kernelDensitySmoothedHistogram(vals, minPlotVal, maxPlotVal, 'color', [1 .7 0]);

% Formatting
set(gcf, 'Position', [34   674   940   206], 'Color', 'w');
title('Principal FC gradient, positive activations')
xlabel('<- Unimodal                Principal FC gradient                Transmodal ->');
hold on
hh = plot_vertical_line(0); set(hh, 'LineStyle', '--');

% add surface
set(gca, 'Position', [0.2700    0.2513    0.7    0.6084])
surfaxhan = axes('Position', [0.0100    0.213    0.2    0.7084]);
surface(pringrad1, 'hcp inflated left', 'nolegend');

drawnow; snapnow;

% Plot negative test_map voxels on same axis
if ~isempty(vals_negative)

    ax_han = gca;
    [density_vals, xValues] = kernelDensitySmoothedHistogram(vals_negative, minPlotVal, maxPlotVal, 'color', [.3 .2 .8]);

    % Formatting
    set(gcf, 'Position', [34   424   940   206], 'Color', 'w');
    title('Principal FC gradient, negative activations')
    xlabel('<- Unimodal                Principal FC gradient                Transmodal ->');
    hold on
    hh = plot_vertical_line(0); set(hh, 'LineStyle', '--');

    % add surface
    set(gca, 'Position', [0.2700    0.2513    0.7    0.6084])
    surfaxhan = axes('Position', [0.0100    0.213    0.2    0.7084]);
    surface(pringrad1, 'hcp inflated left', 'nolegend');

end

drawnow; snapnow;

RESULTS.transuni_vals = [vals, vals_negative];


%% Allen brain project transcriptomic gradients

disp('Allen brain project transcriptomic gradients')

[transcriptomic_grads, transcriptomic_names] = load_image_set('transcriptomic_gradients', 'noverbose');

[vals, vals_negative, minPlotVal, maxPlotVal] = prep_and_extract_values_continuous_reference_maps(transcriptomic_grads, test_map);

axis_handles = plotMatrixKernelDensity(vals, minPlotVal, maxPlotVal, seaborn_colors(), transcriptomic_names);
xlabeltext = 'Transcriptomic gradient value (z-score)';
set(gcf, 'Tag', 'transcriptomics')

% add surfaces and vertical lines at 0
for i = 1:length(transcriptomic_names)

    axes(axis_handles(i))
    hold on
    hh = plot_vertical_line(0); set(hh, 'LineStyle', '--');

    % Move kernel plot axis over
    pos = get(axis_handles(i), 'Position');
    newpos = [pos(1) + .17 pos(2) pos(3) - .12 pos(4)];
    set(axis_handles(i), 'Position', newpos);

    % Create new axis for surface
    surfaxhan = axes('Position', [0.0100    pos(2)    0.2    pos(4)]);
    surface(get_wh_image(transcriptomic_grads, i), 'hcp inflated left', 'nolegend');

end

drawnow; snapnow;

RESULTS.transcriptomic_vals = [vals, vals_negative];


%% Hansen Neuromaps PET neurochemical tracer maps

disp('Neuromaps PET neurochemical tracer maps')

[neurochem, neurochemnames] = load_image_set('hansen22', 'noverbose');

[vals, vals_negative, minPlotVal, maxPlotVal] = prep_and_extract_values_continuous_reference_maps(neurochem, test_map);

axis_handles = plotMatrixKernelDensity(vals, minPlotVal, maxPlotVal, seaborn_colors(size(vals, 2)), neurochemnames);
xlabel('Neurochemical binding value (z-score)');
set(gcf, 'Tag', 'neurotransmitters')

drawnow; snapnow;

% 2nd figure -------------------------------------
create_figure('neurochem', 1, 2);
stats = image_similarity_plot(test_mapz, 'mapset', neurochem, 'networknames', neurochemnames, 'plotstyle', 'polar'); % , 'bicolor');

subplot(1, 2, 2);
plot_colored_bars(stats.r, stats.networknames, seaborn_colors(size(vals, 2)));
set(gca, 'Position', [0.6162    0.110    0.2888    0.8150])
set(gca, 'YDir', 'reverse')

RESULTS.hansen_vals = [vals, vals_negative];


%% Neurosynth topic and term maps

disp('Neurosynth topic and term maps')

ns = load(which('neurosynth_data_obj.mat'));

%        'neurosynth_topics': 
%                   54 topic maps from Yarkoni & Poldrack 2014 topic modeling analysis
%                   Selected from 100 topics for psychological relevance
%                   and given ChatGPT-based summary topic labels by Ke et al. 2024, Nat Neurosci

disp('Neurosynth topics - forward inference maps')
[image_by_feature_correlations, top_feature_tables1, top_ns_maps, bottom_ns_maps, all_ns_r_table] = neurosynth_feature_labels(test_mapz, 'topics_fi', 'images_are_replicates', false, 'noverbose');

disp('Neurosynth topics - reverse inference maps')
[image_by_feature_correlations, top_feature_tables2, top_ns_maps, bottom_ns_maps, all_ns_r_table] = neurosynth_feature_labels(test_mapz, 'topics_ri', 'images_are_replicates', false, 'noverbose');

disp('Neurosynth terms- reverse inference maps')
[image_by_feature_correlations, top_feature_tables3] = neurosynth_feature_labels(test_mapz, 'images_are_replicates', false, 'noverbose');

RESULTS.neurosynth_vals = {top_feature_tables1, top_feature_tables2, top_feature_tables3};

%% Yeo/Buckner lab 

disp('Correlations with Yeo/Bucker resting-state networks')

create_figure('Yeo-Buckner networks', 1, 2);
stats = image_similarity_plot(test_mapz, 'bucknerlab_wholebrain', 'bicolor');

subplot(1, 2, 2)
plot_colored_bars(stats.r, stats.networknames);
set(gca, 'Position', [0.6162    0.300    0.2888    0.5150])

drawnow; snapnow;

% Yeo/Schafer 
[hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(test_mapz, 'atlases', {'yeo17networks'}, 'colors', seaborn_colors(34));

drawnow; snapnow;

RESULTS.networks = output_values_by_region;

end  % main function


%% Subfunctions
% --------------------------------------------


function [vals, vals_negative, minPlotVal, maxPlotVal] = prep_and_extract_values_continuous_reference_maps(reference_map, test_map)
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
reference_mapz = rescale(reference_map, 'zscoreimages');

% make sure voxel lists line up and we haven't removed empty voxels
reference_mapz = replace_empty(reference_mapz);
test_map = replace_empty(test_map);

% binarize and separate into positive and negative values

test_map_binary_pos = test_map;
test_map_binary_pos.dat = single(test_map_binary_pos.dat > 0); % single is conventional format here, but make logical later

test_map_binary_neg = test_map;
test_map_binary_neg.dat = single(test_map_binary_neg.dat < 0); % single is conventional format here, but make logical later

% map the gradient to the test map
% (rather than the reverse, as gradient map is typically less precise
% anyway and interpolation errors matter less)

reference_mapz = resample_space(reference_mapz, test_map_binary_pos);

% Exclude voxels that are exactly 0 or missing in reference map
% This is not appropriate for binary maps
% note: We exclude values that are 0 in all images. If reference maps in
% the set have different missing values in different images (e.g., 0-valued
% maps), this will assume the 0s are valid values unless missing in all
% images in the reference set.)
wh_empty = all(reference_mapz.dat == 0 | isnan(reference_mapz.dat), 2);

% Extract distribution of gradient values for areas in test map, excluding
% missing values in reference map
vals = reference_mapz.dat(logical(test_map_binary_pos.dat & ~wh_empty), :);

vals_negative = reference_mapz.dat(logical(test_map_binary_neg.dat & ~wh_empty), :);

% set limits to plot the histogram
% Make symmetric, and avoid extreme values

% minPlotVal = min(reference_mapz.dat(:));
% maxPlotVal = max(reference_mapz.dat(:));

prcvals = [prctile(reference_mapz.dat(:), 0.1) prctile(reference_mapz.dat(:), 99.9)];

minPlotVal = -max(abs(prcvals));
maxPlotVal = max(abs(prcvals));

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

