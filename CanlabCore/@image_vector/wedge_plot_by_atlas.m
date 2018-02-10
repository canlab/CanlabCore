function hh = wedge_plot_by_atlas(obj_to_plot, varargin)
% Plot a data object or 'signature' pattern divided into local regions
% based on atlas objects.
%
% hh = wedge_plot_by_atlas(obj, varargin)
%
%
% 'signature' : signature mode; changes the statistic that is plotted
% 'atlases' : followed by atlas names as in a cell array defined in load_atlas.m
% 'colors'  : followed by 2-element cell with colors for negative and positive values, respectively
% 'montage' : create montage figure for each atlas
% 'colorband_colors' : followed by cell array (length = number of atlases) with arrays of colors for each region
%
% 'Signature mode':
% Size of wedges: Pattern energy, root mean squared weights in each region, a measure of
% overall importance/contribution to the overall pattern response (assuming that the data the model is applied to is uniform)
%
% Color of wedges: Proportional to valence, cosine similarity to unit
% vector. A measure of whether weights are uniformly positive (red),
% uniformly negative (blue) or mixed (purple).
%
%
% 'Data mode':
% Size of wedges: absolute value of the pattern mean within a parcel
% Color of wedges: sign of the pattern mean within a parcel
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
%  hh = wedge_plot_by_atlas(sig_to_plot,'signature','colors',{[0 0 1],[1 0 0]},'montage')
%
%
% Try with VPS:
% [sigs, signames] = load_image_set('npsplus');
% sig_to_plot = get_wh_image(sigs, 7); signames{7}
% hh = wedge_plot_by_atlas(sig_to_plot, 'signature')
%  hh = wedge_plot_by_atlas(sig_to_plot,'signature','colors',{[0 0 1],[1 0 0]},'montage')


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
%% defaults

%atlases
if ~any(strcmp(varargin,'atlases'))
    atlases = {'basal_ganglia' 'thalamus' }; % 'cerebellum'   'cortex'
else
    atlases = varargin{find(strcmp(varargin,'atlases'))+1};
end
k = length(atlases);

%colors
if ~any(strcmp(varargin,'colors'))
    startcolor = [0 .2 1];
    endcolor = [1 .2 0];
else
    colors=varargin{find(strcmp(varargin,'colors'))+1};
    startcolor=colors{1};
    endcolor=colors{2};
end

% if any(strcmp(varargin,'signature'))
value_limits = [-1 1];
% end

%% initialize cell arrays
mean_weights = {};
sum_sq_weights = {};
signed_msq_weight = {};
msq_weight = {};
vol_in_cubic_mm = {};
voxel_count = {};
atlas_obj = {}; % multiple
labels = {};



%% Compute statistics for plotting

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
    
    [vol_in_cubic_mm{i}, voxel_count{i}] = get_region_volumes(atlas_obj{i});
    
    % divide msq by volume, to get importance per unit volume, add constant
    % to regularize for small regions
    signed_msq_weight{i} = signed_msq_weight{i} ./ (vol_in_cubic_mm{i} + 1000);
    msq_weight{i} = msq_weight{i} ./ (vol_in_cubic_mm{i} + 1000);
    
end

fprintf('Done\n');

%% get colorban colors for wedge plot
% 
if ~any(strcmp(varargin,'colorband_colors'))
    
    for i=1:length(atlases)
        colorband_colors{i} = scn_standard_colors(num_regions(atlas_obj{i}));
    end
    
else
    colorband_colors = varargin{find(strcmp(varargin,'colorband_colors'))+1};

end




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
    
    if any(strcmp(varargin,'signature'))
        data_to_plot = double(msq_weight{i});
        myvalues = pattern_valence{i}'; % get colors
        
    else
        data_to_plot = abs(double(mean_weights{i}));
        myvalues = sign(mean_weights{i}); % get colors
        
    end
    mycolors = values_to_colors(myvalues, value_limits, startcolor, endcolor);
    
    % to-do:  colorband_colors is optional input. pass out colorband_colors
    hh{i} = tor_wedge_plot(data_to_plot, labels{i}, 'colors', mycolors, 'outer_circle_radius', max(abs(data_to_plot)), 'nofigure','colorband','labelstyle','curvy','colorband_colors',colorband_colors{i}); %'bicolor', 'colors', {[1 1 0] [.7 .3 1]},
    
    drawnow
    
end

if any(strcmp(varargin,'signature'))
    disp('Note: Importance is defined as mean squared weight per cubic mm of tissue, regularized by adding a constant (1 cm^3).');
    disp('This is not absolute importance across all voxels (which would favor large regions)')
    disp('And favors small regions with large weights over large regions with local ''hot-spots''.');
    disp('So importance is regularized by adding 1 cm^3 to all region volumes.');
end

%% display regions in separate figure

if any(strcmp(varargin,'montage'))
    
    for i = 1:k
        %fh=create_figure(['Montage of parcellation: ' atlas_obj{i}.atlas_name],1, 1);
        
        o2 = montage(atlas_obj{i}, 'compact2', 'nosymmetric', 'colors',  colorband_colors{i});
% 
%         for r=1:length(o2.activation_maps)
%         colorband_colors{i}{r}=o2.activation_maps{r}.color(1:3);
%         end
        
    end
end



%% todo add legends