function hist_han = histogram(obj, varargin)
% Create a histogram of image values or a series of histograms for each
% image in an image_vector (e.g., fmri_data) object
%
% :Usage:
% ::
%
%     hist_handles = histogram(obj, [optional inputs])
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
%   **obj:**
%       an image_vector object
%
% :Optional Inputs:
%
% **'byimage', 'separate'**
%       Plot histogram for each image separately.
%       The default is to aggregate across all images in the object.
%
% **'by_tissue_type'**
%       Separate into gray, white, CSF tissue compartments and plot
%       histograms for each compartment.
%       This uses the extract_gray_white_csf method, which in turn
%       currently uses the images
%       'gray_matter_mask.img' 'canonical_white_matter.img' 'canonical_ventricles.img'
%       These images are based on the SPM8 a priori tissue probability
%       maps, but they have been cleaned up and made symmetrical and/or eroded
%       so that the white and CSF compartments are unlikely to contain very
%       much gray matter.  The gray compartment is currently more
%       inclusive. The potential value of this is that signal in the CSF/white
%       compartments may be removed from images prior to/during analysis
%
% **'plot'**
%       No function right now.
%
% **'nbins'**
%       Followed by number of bins to use (default 100)
%
% **'color'**
%       Followed by color [r g b] triplet or text, e.g., 'b'
%
% **'mask'**
%       Followed by mask image name or fmri_data object
%
% **'doline'**
%       Plot lines for histograms
%
% **'noline'**
%       Suppress lines for histograms
%
% **'nofill'**
%       Suppress histogram fill
%
% **'nofigure'**
%       Suppress figure creation
%
% :Outputs:
%
%   **hist_han:**
%        Graphics handles to histogram lines or bars
%
%   **out2:**
%        description of out2
%
% :Examples:
% ::
%
% hist_han = histogram(wp_alone, 'byimage', 'color', 'b');
% hist_han = histogram(wp_alone, 'byimage', 'by_tissue_type');
% hist_han = histogram(wp_alone, 'by_tissue_type');
%
%
% :References:
%   None. These are just histograms.
%
% :See also:
%   - image_histogram, image_intensity_histograms
%
%
% ..
%    Programmers' notes:
%   Documentation and update, July 2016, Tor Wager
% ..

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

doplot = 1;
doline = 1;
dofill = 1;
dofigure = 1;
do_by_image = 0;
nbins = [];
color = [.3 .3 .3];
mask = [];
by_tissue_type = 0;

% initalize optional variables to default values here.


% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'byimage', 'by_image', 'separate'}, do_by_image = 1;
                
            case 'plot', doplot = 1;
                
            case 'nbins', nbins = varargin{i+1}; varargin{i+1} = [];
            case 'color', color = varargin{i+1}; varargin{i+1} = [];
            case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
                
            case 'doline', doline = 1;
            case 'noline', doline = 0;
            case 'nofill', dofill = 0;
            case {'nofigure', 'nofig'}, dofigure = 0;
                
            case 'by_tissue_type', by_tissue_type = 1; varargin{i} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% -------------------------------------------------------------------------
% SPECIAL FUNCTION: BY TISSUE TYPE
% -------------------------------------------------------------------------
if by_tissue_type
    
    [values, components, full_data_objects] = extract_gray_white_csf(obj);
    
    hist_han = cell(1, 3);
    colors = {[.2 .2 .2] [1 .6 .6] [.3 .3 1]};
    
    for i = 1:3
        
        hist_han{i} = histogram(full_data_objects{i}, varargin{:}, 'color', colors{i}, 'nofigure');
        
        full_data_objects{i} = remove_empty(full_data_objects{i});  % just in case
        
        stdevs(:, i) = nanstd(full_data_objects{i}.dat);
        
    end
        
    
    % Plot locations of each and relationships
     
    create_figure('relationships', 1, 3);
    disp('------------------------------------------------------------')
    disp('Plotting mean values by tissue type')
    fprintf('Gray = gray, red = white, blue = CSF\n');

    barplot_columns(values, 'nofig', 'colors', colors, 'names', {'gray' 'white' 'CSF'});
    %set(gca, 'XTickLabel', {'gray' 'white' 'CSF'});
    title('Mean values by tissue type');
    
    subplot(1, 3, 2);
    disp('------------------------------------------------------------')
    disp('Plotting spatial standard deviations across voxels by tissue type')
    fprintf('Gray = gray, red = white, blue = CSF\n');
    
    barplot_columns(stdevs, 'nofig', 'colors', colors, 'names', {'gray' 'white' 'CSF'});
    %set(gca, 'XTickLabel', {'gray' 'white' 'CSF'});
    title('Std. Dev. by tissue type');
    
    subplot(1, 3, 3);
    plotmatrix([values stdevs]);
    title('Relations [means stdevs]')
    
    return
    
end

% -------------------------------------------------------------------------
% MASKING
% -------------------------------------------------------------------------

if ~isempty(mask)
    
    if isa(mask, 'image_vector')
        obj = apply_mask(obj, mask);
    else
        obj = apply_mask(obj, fmri_data(mask));
    end
    
end

% -------------------------------------------------------------------------
% HISTOGRAMS
% -------------------------------------------------------------------------

if do_by_image
    % By image: one figure per image in object
    % -----------------------------------------------------
    
    % Set scale range
    Xtmp = obj.dat(:);
    Xtmp(Xtmp == 0 | isnan(Xtmp)) = [];
    mylim = prctile(Xtmp, [1 99]);
    
    % Set n bins if auto-select
    if isempty(nbins), nbins = ceil(size(obj.dat, 1) ./ 500); end
    
    % Set up subplots
    nimgs = size(obj.dat, 2);
    [rows, cols] = deal(floor(nimgs .^ .5));
    while rows * cols < nimgs, cols = cols + 1; end
    naxes = rows * cols;
    
    if dofigure
        create_figure('histogram', rows, cols);
    end
    
    % Do for each image
    for i = 1:nimgs
        
        subplot(rows, cols, i);
        hist_han(i) = create_hist(obj.dat(:, i), nbins, doline, dofill, color);
        
        % set axes
        set(gca, 'XLim', mylim);
        
    end
    
    for i = nimgs + 1:naxes
        subplot(rows, cols, i)
        axis off
    end
    
else
    % Overall: one histogram
    % -----------------------------------------------------
    
    % Set n bins if auto-select
    if isempty(nbins), nbins = ceil(size(obj.dat, 1) ./ 500); end
    
    if dofigure
        create_figure('histogram');
    end
    
    Xtmp = obj.dat(:);
    hist_han = create_hist(Xtmp, nbins, doline, dofill, color);
    
    % Line at 0
    hh = plot_vertical_line(0);
    set(hh, 'LineStyle', '--');

    xlabel('Values'); ylabel('Frequency');
    title('Histogram of values');
    
end





end % FUNCTION


function han = create_hist(Xtmp, nbins, doline, dofill, color)

% remove zeros
Xtmp(Xtmp == 0 | isnan(Xtmp)) = [];

% Get bins
x = linspace(prctile(Xtmp, .01), prctile(Xtmp, 99.9), nbins);

[h, x] = hist(Xtmp, x);

% convert to PDF
h = h ./ sum(h);

if doline
    han = plot(x, h, 'LineWidth', 3, 'color', color);
    
else
    han = bar(x, h);
    set(han, 'FaceColor', color, 'EdgeColor', 'none');
    
end

if dofill
    fill([x(1)-eps x x(end)+eps], [0 h 0], color, 'FaceAlpha', .2);
end

hold on; % still trouble with automatic turn-off
 
axis tight

% set x-axis
% should not need to do this if we've chosen bins carefully above
% wh = find(h > .0001);
% xlims = [x(min(wh)) x(max(wh))];
% set(gca, 'XLim', xlims);

drawnow

end


