function riverplot(layer1fmri_obj, varargin)
% Make a riverplot of relationships among images, stored in an fmri_data object,
% or two sets of images, stored in two fmri_data objects.
%
% - With one input, uses self
% - With two inputs, uses first as Layer 1 and second as Layer 2 of riverplot
% - Uses obj.image_names for titles of images, so set these before calling this function
%
% :Usage:
% ::
%
%     [list outputs here] = function_name(list inputs here, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
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
%   **layer1fmri_obj:**
%        An fmri_data object.  Default is to plot associations between
%        images in this object.
%
% :Optional Inputs:
%   **layer2:**
%        A second fmri_data object.  This is layer 2 of the plot.  Plots
%        assocations between images in layer1 object and layer2 object
%
%   **r:**
%        Change similarity metric to correlation.  Default is cosine
%        similarity.
%
%   **dice:**
%        Change similarity metric to dice coefficient.  Default is cosine
%        similarity.
%
%   **pos:**
%        Change similarity matrix thresholding to positive values only.
%        Default is to use original similarity matrix.
%
%   **neg:**
%        Change similarity matrix thresholding to negative values only.
%        Default is to use original similarity matrix.
%
%   **abs:**
%        Change similarity matrix thresholding to absolute values.
%        Default is to use original similarity matrix.
%
%   **reorder:**
%        Automatic reordering of rows of plot to try to improve visibility.
%
%   **recolor:**
%        Automatic color assignment of layer2 colors based on layer1
%        weights
%
%   **layer1order:**
%        Followed by custom order for Layer 1, specified by vector of
%        integers, e.g., [1 3 2 4 5 6]
%
%   **layer2order:**
%        Followed by custom order for Layer 2, specified by vector of
%        integers, e.g., [1 3 2 4 5 6]
%
%   **colors1:**
%        Followed by cell vector of [r g b] colors for each image in Layer 1
%
%   **colors2:**
%        Followed by cell vector of [r g b] colors for each image in Layer 2
%
%   **thin:**
%        Turn off thick/full-rectangle ribbon width coming from Layer 1 
%
%
% :Outputs:
%
%   **layer1:**
%        Vector of structures describing rectangles in Layer 1, with
%        handles, positions, etc.
%
%   **layer2:**
%        Vector of structures describing rectangles in Layer 2, with
%        handles, positions, etc.
%
%   **ribbons:**
%        Vector of structures describing ribbons, with
%        handles, positions, etc.
%
%   **layer2colors:**
%        cell array specifying layer2 colors if recolor option is used
%        
%
% :Examples:
% ::
%
% layer1fmri_obj = fmri_data(placeboimgs);
% layer2fmri_obj = fmri_data(appraisalimgs);
% riverplot(layer1fmri_obj, 'layer2', layer2fmri_obj);
%
% Add colors and names:
% layer1colors = {[.8 .8 .3] [.3 .3 1]};
% layer2colors = {[0 .2 .8] [1 1 .2] [1 .4 .4] [0 .6 0] [1 .7 0]};  %seaborn_colors(n1);
%
% layer1names = {'Placebo increases', 'Placebo decreases'};  %stats.networknames;
% layer2names = {'Neg Emo' 'Pos Emo' 'Emo Reg' 'Social Cog' 'Self'};
%
% layer1fmri_obj.image_names = char(layer1names{:});
% layer2fmri_obj.image_names = char(layer2names{:});
%
% riverplot(layer1fmri_obj, 'layer2', layer2fmri_obj, 'layer1colors', layer1colors, 'layer2colors', layer2colors);
%
% Example: Plot Buckner lab network maps with themselves:
% [bucknermaps, networknames] = load_image_set('bucknerlab');
% bucknermaps.image_names = char(networknames{:});
% 
% layer1colors = bucknerlab_colors(7);
% riverplot(bucknermaps, 'pos', 'reorder', 'layer1colors', layer1colors, 'layer2colors', layer1colors);
% 
% riverplot(bucknermaps, 'pos', 'layer1colors', layer1colors, 'layer2colors', layer1colors);
% 
%
% :References:
%   Based off of concept in River Plot R package.
%
% :See also:
%   - image_similarity_plot, tor_polar_plot
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%    Created July 2016 by Tor Wager
% ..

% BELOW IS A STANDARD TEMPLATE FOR DEFINING VARIABLE (OPTIONAL) INPUT
% ARGUMENTS. MANY FUNCTIONS NEED TO PARSE OPTIONAL ARGS, SO THIS MAY BE
% USEFUL.


% ..
%    DEFAULTS AND INPUTS
% ..



% optional inputs with default values

layer2fmri_obj = layer1fmri_obj;  % Default: Plot object with itself
sim_metric = 'cosine_sim';
threshold = 'none';

nlayer1 = size(layer1fmri_obj.dat, 2);
nlayer2 = size(layer2fmri_obj.dat, 2);

layer1colors = custom_colors([.2 .8 .2], [.8 .8 .2], nlayer1);
layer2colors = custom_colors([.2 .2 .8], [.2 .8 .8], nlayer2);

doreorder = 0;
recolor = 0;

layer1coveragestr = 'layer1fullcoverage';

layer1order = [];

layer2order = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'layer2', 'layer2fmri_obj'},
                layer2fmri_obj = varargin{i+1}; varargin{i+1} = [];
                nlayer2 = size(layer2fmri_obj.dat, 2);
                layer2colors = custom_colors([.2 .2 .8], [.2 .8 .8], nlayer2);
                
            case 'sim_metric', sim_metric = varargin{i+1}; varargin{i+1} = [];
            case 'threshold', threshold = varargin{i+1}; varargin{i+1} = [];
                
            case {'colors1', 'layer1colors'}, layer1colors = varargin{i+1}; varargin{i+1} = [];
            case {'colors2', 'layer2colors'}, layer2colors = varargin{i+1}; varargin{i+1} = [];
                
            case 'reorder', doreorder = 1;
              case 'recolor', recolor = 1;
                  
            case {'r', 'corr', 'correlation'}
                sim_metric = 'r';
                
            case 'dice'
                sim_metric = 'dice';
                
            case 'pos'
                threshold = 'pos';
                
            case 'neg'
                threshold = 'neg';
                
            case 'abs'
                threshold = 'abs';
                
            case 'thin'
                layer1coveragestr = ' ';

            case 'layer1order', layer1order = varargin{i+1}; varargin{i+1} = [];
            case 'layer2order', layer2order = varargin{i+1}; varargin{i+1} = [];

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% Get similarity matrix
% -------------------------------------------------------------------------

switch sim_metric
    
    case 'cosine_sim' % default
        
        stats = image_similarity_plot(layer1fmri_obj, 'mapset', layer2fmri_obj, 'noplot', 'cosine_similarity');
        sim_matrix = stats.r;
        
    case {'r', 'corr', 'correlation'}
        
        stats = image_similarity_plot(layer1fmri_obj, 'mapset', layer2fmri_obj, 'noplot');
        sim_matrix = stats.r;
        
    case 'dice'
        
        dice_coeff = dice_coeff_image(cat(layer1fmri_obj, layer2fmri_obj));
        
        sim_matrix = dice_coeff(1:nlayer1, nlayer1+1:end)';
        
    otherwise
        error('Unknown sim_metric. Choose r, cosine_sim, or dice');
end

% -------------------------------------------------------------------------
% Threshold
% -------------------------------------------------------------------------

switch threshold
    
    case 'none'
        % do nothing - default
        
    case {'pos', 'positive'}
        
        % Threshold to look at positive associations
        sim_matrix(sim_matrix < 0) = 0;
        
    case {'neg', 'negative'}
        sim_matrix = -sim_matrix;
        sim_matrix(sim_matrix < 0) = 0;
        
    case 'abs'
        sim_matrix = abs(sim_matrix);
        
    otherwise
        error('Unknown threshold type. Choose pos, neg, or abs');
end

% -------------------------------------------------------------------------
% Names and sizes
% -------------------------------------------------------------------------

layer1names = format_strings_for_legend(layer1fmri_obj.image_names);
layer2names = format_strings_for_legend(layer2fmri_obj.image_names);

[n2, n1] = size(sim_matrix);

if length(layer1names) < n1, error('Layer 1: Not enough names entered. Add names in input layer1 fmri_object.image_names'); end
if length(layer2names) < n2, error('Layer 2: Not enough names entered. Add names in input layer2 fmri_object.image_names'); end

% -------------------------------------------------------------------------
% Reorder rows and columns if asked for
% -------------------------------------------------------------------------

if doreorder
    
    %  REORDER rows here based on similarity to make plot look cleaner
    [sim_matrix, layer2names, layer2colors] = riverplot_reorder_matrix(sim_matrix, layer2names, layer2colors);
    
end

if ~isempty(layer1order)
    
    sim_matrix = sim_matrix(:, layer1order);
    layer1names = layer1names(layer1order);
    layer1colors = layer1colors(layer1order);
    
end

if ~isempty(layer2order)
    
    sim_matrix = sim_matrix(layer2order, :);
    layer2names = layer2names(layer2order);
    layer2colors = layer2colors(layer2order);
    
end

if recolor
    layer2colors=riverplot_recolor_layer2(sim_matrix,layer1colors);
end

% -------------------------------------------------------------------------
% Find start coordinate points for layers
% -------------------------------------------------------------------------

layer2x = max(n1, n2) ./ 2;  % 1/2 largest array size, to keep x proportions balanced with y

% each rect is 1 unit, starting at y = 0 for layer. so if same size, no y_loc offset
% midpoint should be same for both
mid2 = (n2 + 1) ./ 2;
% mid1 should equal mid2  : % (n1 + 1) ./ 2 - y_loc = mid 2
y_loc = -(n1 + 1) ./ 2 + mid2;

% -------------------------------------------------------------------------
% Create the plot
% -------------------------------------------------------------------------

create_figure('riverplot');
set(gca, 'YDir', 'reverse');
layer1 = riverplot_draw_layer(0, n1, 'colors', layer1colors, 'y_loc', y_loc);
layer2 = riverplot_draw_layer(layer2x, n2, 'colors', layer2colors);

ribbons = riverplot_ribbon_matrix(layer1, layer2, sim_matrix, 'colors', layer1colors, layer1coveragestr, 'steepness', 0);

% Turn off lines
riverplot_toggle_lines(ribbons);

% Increase opacity
riverplot_set_ribbon_property(ribbons, 'FaceAlpha', .6);

set(gca, 'XLim', [-3 layer2x+3]);

% Add names
layer1 = riverplot_layer_names(layer1, layer1names);
layer2 = riverplot_layer_names(layer2, layer2names, 'right');

axis off

end % function
