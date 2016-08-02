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
%       Plot histogram for each image separately
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
nbins = 100;
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
    
    fprintf('Gray = gray, red = white, blue = CSF\n');
    
    
    % Plot locations of each and relationships
     
    create_figure('relationships', 1, 3);
    barplot_columns(values, 'nofig', 'colors', colors);
    set(gca, 'XTickLabel', {'gray' 'white' 'CSF'});
    title('Means');
    
    subplot(1, 3, 2);
    barplot_columns(stdevs, 'nofig', 'colors', colors);
    set(gca, 'XTickLabel', {'gray' 'white' 'CSF'});
    title('Std. Dev.');
    
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
    
    if dofigure
        create_figure('histogram');
    end
    
    Xtmp = obj.dat(:);
    hist_han = create_hist(Xtmp, nbins, doline, dofill, color);
    
    xlabel('Values'); ylabel('Frequency');
    title('Histogram of values');
    
end





end % FUNCTION


function han = create_hist(Xtmp, nbins, doline, dofill, color)

% remove zeros
Xtmp(Xtmp == 0 | isnan(Xtmp)) = [];

[h, x] = hist(Xtmp, nbins);

% convert to PDF
h = h ./ sum(h);

if doline
    han = plot(x, h, 'LineWidth', 3, 'color', color);
    
else
    han = bar(x, h);
    set(han, 'FaceColor', color, 'EdgeColor', 'none');
    
end

if dofill
    fill(x, h, color, 'FaceAlpha', .2);
end

hold on; % still trouble with automatic turn-off
 
axis tight

drawnow

end


