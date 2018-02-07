function hh = wedge_plot_by_atlas(obj_to_plot, varargin)
% Plot a data object or 'signature' pattern divided into local regions
% based on atlas objects.
%
% hh = wedge_plot_by_atlas(obj, varargin)
%
%
% 'signature' : signature mode; change what is plotted
% 'atlases' : followed by atlas names as defined in load_atlas.m
% 'colors'  : followed by 2-element cell with colors for negative and positive values, respectively
%
% 'Signature mode': 
% Size of wedges: Pattern energy, root mean squared weights in each region, a measure of
% overall importance/contribution to the overall pattern response
%
% Color of wedges: Proportional to valence, cosine similarity to unit
% vector. A measure of whether weights are uniformly positive (red),
% uniformly negative (blue) or mixed (purple).
% 
% Examples:
% 
% Load sig
% nps = load_image_set('npsplus');
% nps = get_wh_image(nps, 1);
% sig_to_plot = nps;
%
% Plot:
%  hh = wedge_plot_by_atlas(sig_to_plot, 'signature')
%
% Try with VPS:
% [sigs, signames] = load_image_set('npsplus');
% sig_to_plot = get_wh_image(sigs, 7); signames{7}
% hh = wedge_plot_by_atlas(sig_to_plot, 'signature')


% Notes on pattern valence
% ------------------------------------------------------------------
% built into apply_parcellation
%
% Define pattern valence as similarity with unit vector (all positive,
% identical weights)
%
% If x is a vector of pattern weights, 
% The cosine similarity with the unit vector is a measure of how uniform
% the weights are, on a scale of 1 to -1.  1 indicates that the pattern
% computes the region average (all weights identical and positive). -1
% indicates that the pattern computes the negative region average (all
% weights identical and negative).  



%
%%
atlases = {'basal_ganglia' 'thalamus' 'cerebellum'}; % 'cortex' 

k = length(atlases);

mean_weights = {};
sum_sq_weights = {};
signed_msq_weight = {};
msq_weight = {};
vol_in_cubic_mm = {};
voxel_count = {};
atlas_obj = {}; % multiple
labels = {};

% Prep
disp('Prepping atlases: ');

for i = 1:k
        
    fprintf('%s ', atlases{i});
    
    atlas_obj{i} = load_atlas(atlases{i});
    
    % redefine names if needed
    switch atlases{i}
        case 'cerebellum'
            atlas_obj{i}.labels = strrep(atlas_obj{i}.labels, 'Left', 'L');
            atlas_obj{i}.labels = strrep(atlas_obj{i}.labels, 'Right', 'R');
            atlas_obj{i}.labels = strrep(atlas_obj{i}.labels, 'Vermis', 'Verm');
            atlas_obj{i}.labels = strrep(atlas_obj{i}.labels, 'Interposed', 'Intp');
            atlas_obj{i}.labels = strrep(atlas_obj{i}.labels, 'Dentate', 'Dnt');
            atlas_obj{i}.labels = strrep(atlas_obj{i}.labels, 'Fastigial', 'Fst');
    end
    
    labels{i} = format_strings_for_legend(atlas_obj{i}.labels);
    
    [mean_weights{i}, sum_sq_weights{i}, pattern_valence{i}] = apply_parcellation(obj_to_plot, atlas_obj{i}, 'pattern_expression', obj_to_plot);
    
    signed_msq_weight{i} = sign(mean_weights{i}) .* sum_sq_weights{i} .^ .5;

    msq_weight{i} = sum_sq_weights{i} .^ .5;
    
    [vol_in_cubic_mm{i} voxel_count{i}] = get_region_volumes(atlas_obj{i});
    
    % divide msq by volume, to get importance per unit volume, add constant
    % to regularize for small regions
    signed_msq_weight{i} = signed_msq_weight{i} ./ (vol_in_cubic_mm{i} + 1000);
    msq_weight{i} = msq_weight{i} ./ (vol_in_cubic_mm{i} + 1000);
    
end

fprintf('Done\n');

%% Get colors based on pattern valence. Define consistent elemens across atlases:

value_limits = [-1 1];
startcolor = [0 .2 1]; 
endcolor = [1 .2 0];

%% Plot importance (squared weights, regularized), 
% with colors indicating valence (whether weights are uniformly pos (red) or neg (blue)
% mixed weights: purple.

disp('Note: Mean weights reflect homogeneity in sign and magnitude across region,')
disp('not high spatial frequency/pattern information.');


create_figure('wedge overall importance', 1, k)
% overall importance of region is proportional to geometric mean of weights

clear hh

for i = 1:k
    
    subplot(1, k, i);
 
    data_to_plot = double(msq_weight{i});                
    
        % get colors
    myvalues = pattern_valence{i}';
    mycolors = values_to_colors(myvalues, value_limits, startcolor, endcolor);


    hh{i} = tor_wedge_plot(data_to_plot, labels{i}, 'colors', mycolors, 'outer_circle_radius', max(abs(data_to_plot)), 'nofigure'); %'bicolor', 'colors', {[1 1 0] [.7 .3 1]},
    
    drawnow

end

disp('Note: Importance is defined as mean squared weight per cubic mm of tissue, regularized by adding a constant (1 cm^3).');
disp('This is not absolute importance across all voxels (which would favor large regions)')
disp('And favors small regions with large weights over large regions with local ''hot-spots''.');
disp('So importance is regularized by adding 1 cm^3 to all region volumes.');

% to do:
% display regions

