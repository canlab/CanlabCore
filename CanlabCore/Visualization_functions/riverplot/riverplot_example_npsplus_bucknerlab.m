% This example loads the 'nps plus' image set of patterns, and examines the
% relationships between these patterns and the standard Buckner Lab
% networks.

% Layer 1 of the river plot is the Wager lab patterns, and Layer 2 is the
% Buckner Lab networks

%% Load NPSplus

[npsplus_obj, layer1names, imgnames] = load_image_set('npsplus');

wh = 1:4;
npsplus_obj = get_wh_image(npsplus_obj, wh);
layer1names = layer1names(wh);

layer1names

%% Get rectangular similarity matrix

stats = image_similarity_plot(npsplus_obj, 'bucknerlab', 'noplot', 'cosine_similarity');
% Similarity matrix is stats.r

% Threshold to look at positive associations
sim_matrix = stats.r;
sim_matrix(sim_matrix < 0) = 0;

layer2names = stats.networknames;

[n2, n1] = size(sim_matrix);

%layer1colors = scn_standard_colors(n1);
%layer1colors = seaborn_colors(n1);

layer1colors = {[1 .4 .2] [.3 .3 1] [.7 .3 1] [.3 .7 .5]};

layer2colors = seaborn_colors(n2);

%% Reorder rows and columns for display, as needed

% Custom REORDERING of layer 1 for aesthetics
neworder = [1 3 2 4];
sim_matrix = sim_matrix(:, neworder);
layer1names = layer1names(neworder);
layer1colors = layer1colors(neworder);

%  REORDER rows here based on similarity to make plot look cleaner
[sim_matrix, layer2names, layer2colors] = riverplot_reorder_matrix(sim_matrix, layer2names, layer2colors);

% Custom REORDERING of layer 2 because the automatic one isn't super
neworder = [1 7 2 3 4 5 6];
sim_matrix = sim_matrix(neworder, :);
layer2names = layer2names(neworder);
layer2colors = layer2colors(neworder);

%% Create the plot

create_figure('riverplot'); 
set(gca, 'YDir', 'reverse');
layer1 = riverplot_draw_layer(2, n1, 'colors', layer1colors, 'y_loc', 1.5);
layer2 = riverplot_draw_layer(5, n2, 'colors', layer2colors);

ribbons = riverplot_ribbon_matrix(layer1, layer2, sim_matrix, 'colors', layer1colors, 'layer1fullcoverage', 'steepness', 0);

% Turn off lines
riverplot_toggle_lines(ribbons);

% Increase opacity
riverplot_set_ribbon_property(ribbons, 'FaceAlpha', .6);

set(gca, 'XLim', [0 8]);

% Add names
layer1 = riverplot_layer_names(layer1, layer1names);
layer2 = riverplot_layer_names(layer2, layer2names, 'right');

axis off

