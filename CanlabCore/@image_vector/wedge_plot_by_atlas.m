function [hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(obj_to_plot, varargin)
% Plot a data object or 'signature' pattern divided into local regions
% based on atlas objects.
%
% [hh, output_values_by_region, labels] = wedge_plot_by_atlas(obj, varargin)
%
%
% 'signature' : signature mode; changes the statistic that is plotted
% 'atlases' : followed by atlas names as in a cell array defined in load_atlas.m
% 'colors'  : followed by 2-element cell with colors for negative and positive values, respectively
% 'montage' : create montage figure for each atlas
% 'colorband_colors' : followed by cell array (length = number of atlases) with arrays of colors for each region
%
% 'Signature mode':
% Size of wedges: Pattern energy, related to root mean squared weights in each region
%                 Normalized local pattern expression (sqrt(w''*w)/(vol_in_mm+1000)
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
% Outputs:
% ------------------------------------------------------------------------
% output_values_by_region: Cell array, one cell per atlas, with a matrix of
% images x regions in atlas. Values depend on whether you use 'signature
% mode' or 'data mode'.
%
% labels: region names for each atlas
%
% atlas_obj: cell array of atlas objects
%
% colorband_colors: cell array of colors for each atlas
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
%
% Try some data:
% imgs = load_image_set('emotionreg');
% hh = wedge_plot_by_atlas(imgs, 'atlases', {'cit168' 'brainstem'});
%
% Try custom colors, mirroring clusters across L/R hem networks:
% [colors1, colors2] = deal(scn_standard_colors(16));
% colors = {};
% indx = 1;
% for i = 1:length(colors1)
%     colors{indx} = colors1{i}; colors{indx + 1} = colors2{i}; indx = indx + 2;
% end
% [hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(imgs, 'atlases', {'yeo17networks'}, 'montage', 'colorband_colors', colors);

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

% Signature mode
do_sig_mode = false;
if any(strcmp(varargin, 'signature')), do_sig_mode = true; end

% atlases
if ~any(strcmp(varargin,'atlases'))
    atlases = {'basal_ganglia' 'thalamus' }; % 'cerebellum'   'cortex'
else
    atlases = varargin{find(strcmp(varargin,'atlases'))+1};
end
k = length(atlases);

% colors
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
    
    if do_sig_mode
        
        [mean_weights{i}, sum_sq_weights{i}, pattern_valence{i}] = apply_parcellation(obj_to_plot, atlas_obj{i}, 'pattern_expression', obj_to_plot);
        
        signed_msq_weight{i} = sign(mean_weights{i}) .* sum_sq_weights{i} .^ .5;
        
        msq_weight{i} = sum_sq_weights{i} .^ .5;
        
        [vol_in_cubic_mm{i}, voxel_count{i}] = get_region_volumes(atlas_obj{i});
        
        % divide msq by volume, to get importance per unit volume, add constant
        % to regularize for small regions
        signed_msq_weight{i} = signed_msq_weight{i} ./ (vol_in_cubic_mm{i} + 1000);
        msq_weight{i} = msq_weight{i} ./ (vol_in_cubic_mm{i} + 1000);
        
        output_values_by_region = msq_weight;
        
        legendtext = sprintf('Wedge plots depict normalized local pattern expression (sqrt(w''*w)/(vol_in_mm+1000) based on the signature pattern visualized. ');
        legendtext = [legendtext 'The color of each wedge indicates the pattern valence, whether values are positive (red) or negative (blue) on average. '];
        legendtext = [legendtext 'A region may have high pattern importance (large weights, indicated by a large wedge) and either net positive weights (red), negative weights (blue), or mixed weights (purple). If weights are homogenous (bright red or blue), the region average is a good summary of the pattern. But if weights are mixed, the region average is not a good summary of the pattern. '];
        
        legendtext = [legendtext '\nNote: Importance is defined as mean squared weight per cubic mm of tissue, regularized by adding a constant (1 cm^3). This is not absolute importance across all voxels (which would favor large regions and favors small regions with large weights over large regions with local ''hot-spots''. Therefore, importance is regularized by adding 1 cm^3 to all region volumes'];

    else
        % Data mode
        
        mean_weights{i} = apply_parcellation(obj_to_plot, atlas_obj{i});
        
        output_values_by_region = mean_weights;
        
        legendtext = sprintf('Wedge plots depict mean images values across voxels. Red indicates positive values and blue negative values. If multiple images were entered, the darker shaded area indicates the standard error of the mean (SEM) across individuals.');

    end
    
    
end


%% get colorband colors for wedge plot
% 
if ~any(strcmp(varargin,'colorband_colors'))
    
    for i=1:length(atlases)
        colorband_colors{i} = scn_standard_colors(num_regions(atlas_obj{i}));
    end
    
else
    
    for i=1:length(atlases)
        colorband_colors{i} = varargin{find(strcmp(varargin,'colorband_colors'))+1};
    end
    
    if ~iscell(colorband_colors), colorband_colors = {colorband_colors}; end

end



%% Plot importance (squared weights, regularized),
% with colors indicating valence (whether weights are uniformly pos (red) or neg (blue)
% mixed weights: purple.

disp('Note: Mean weights reflect homogeneity in sign and magnitude across region,')
disp('not high spatial frequency/pattern information.');


create_figure('wedge overall importance', 1, k);
% overall importance of region is proportional to geometric mean of weights

clear hh

for i = 1:k % k indexes atlases. so do this for each atlas.
    
    subplot(1, k, i);
    
    if length(colorband_colors{i}) ~= size(mean_weights{i}, 2)
        disp('Number of colorband colors must equal number of regions displayed.')
        disp('Create one cell containing a cell array of colors for each atlas.')
        disp('Pausing so you can check:');
        keyboard
    end
        
    if do_sig_mode
        
        data_to_plot = double(msq_weight{i});
        myvalues = pattern_valence{i}'; % get colors
        myouterradius = max(abs(data_to_plot));
        
        mycolors = values_to_colors(myvalues, value_limits, startcolor, endcolor);
        
        % to-do:  colorband_colors is optional input. pass out colorband_colors
        hh{i} = tor_wedge_plot(data_to_plot, labels{i}, 'colors', mycolors, 'outer_circle_radius', myouterradius, 'nofigure','colorband','labelstyle','curvy','colorband_colors', colorband_colors{i}); %'bicolor', 'colors', {[1 1 0] [.7 .3 1]},
    
    else % data mode
        
        data_to_plot = double(mean_weights{i});   % mean weights really average values, not mean weights
        mymean = mean(data_to_plot, 1);
        if size(data_to_plot,1) > 1
            myste = ste(data_to_plot);
        else
            myste = 0;
        end
        
        myouterradius = max(abs(mymean) + myste);
        mycolors = {[1 0 .2] [.2 0 1]}; % {[1 1 0] [.7 .3 1]}
        
        hh{i} = tor_wedge_plot(data_to_plot, labels{i}, 'outer_circle_radius', myouterradius, 'nofigure','colorband','labelstyle','curvy','colorband_colors', colorband_colors{i}, 'bicolor', 'colors', mycolors);
        
    end

    drawnow
    
end

% Make the figure big
%set(gcf, 'Units', 'Normalized', 'Position', [.05 .3 .5 .6])

% Descriptive legend text - enter a series of cell arrays with text
% --------------------------------------------------------------
headerstr = 'Wedge plot: ';
canlab_print_legend_text(headerstr, legendtext);

%% display regions in separate figure

if any(strcmp(varargin,'montage'))
    
    for i = 1:k
        %fh=create_figure(['Montage of parcellation: ' atlas_obj{i}.atlas_name],1, 1);
        
        % colorband_colors specifies colors for each region. Each cell
        % should contain a cell array with colors  length(num regions in atlas_obj{i})

        o2 = montage(atlas_obj{i}, 'compact2', 'nosymmetric', 'colors',  colorband_colors{i});
        drawnow, snapnow
        
        % Individual regions
        
         montage(atlas_obj{i}, 'nosymmetric', 'regioncenters', 'colors',  colorband_colors{i});
         
         set(gcf, 'Tag', 'Individual regions'); % change tag to avoid figure resizing on repeated calls
         drawnow, snapnow
         
    end
end



%% todo add legends


end % function
