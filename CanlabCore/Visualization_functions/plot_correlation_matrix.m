function OUT = plot_correlation_matrix(X, varargin)
% Plots an image of a correlation matrix with circles, text, image, or combinations
%
% :Usage:
% ::
%
%     OUT = plot_correlation_matrix(X, [optional inputs])
%
% - Can do full or partial correlations, Spearman or Pearson's
% - FDR correction across pairwise tests if desired.
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Tor Wager
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
%   **X:**
%        n_observations x k_variables data matrix
%
% :Optional Inputs:
%   **var_names:**
%        Followed by cell array of variable names
%
%   **p_thr:**
%        Followed by P-value threshold for significant correlations
%
%   **[OTHERS]:**
%       See code for other optional inputs controlling display
%       Enter a keyword followed by a value (e.g., true / false)
%       These include:
%       dofigure, dospearman, doimage, docircles, dotext, colorlimit, text_x_offset, text_y_offset,  
%       text_fsize, text_nonsig_color, text_sig_color
%       'dopartitions', followed by k-length integer vector defining partitions of variables  [optional] 
%       'partitioncolors', m-length cell vector defining rgb color triplets for each partition [optional] 
%       'p_thr', followed by p-value threshold
%
% :Outputs:
%
%   **OUT:**
%        structure with various outputs
%
% :Examples:
% ::
% % Generate some simulated data (N = 50 cases) and plot correlations: 
% S = toeplitz([1 .6 .3 .1 0 0]);
% X = mvnrnd([0 0 0 0 0 0], S, 50);
% var_names = {'A' 'B' 'C' 'D' 'E' 'F'};
% OUT = plot_correlation_matrix(X, 'var_names', var_names);
%
% % Also, e.g.:
%   OUT = plot_correlation_matrix(X, 'var_names', var_names, 'colorlimit', [-.5 .5]);
%
% % For larger matrices, do not plot circles:
% OUT = plot_correlation_matrix(X, 'doimage', true, 'docircles', false);
%
% % Can do full or partial correlations, Spearman or Pearson's
% % FDR correction across pairwise tests if desired.
% OUT = plot_correlation_matrix(X, 'doimage', true, 'docircles', false, 'dospearman', true, 'dopartial', true, 'dofdr', true);
%
% :References:
%   None - basic stats functions.
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% BELOW IS A STANDARD TEMPLATE FOR DEFINING VARIABLE (OPTIONAL) INPUT
% ARGUMENTS. MANY FUNCTIONS NEED TO PARSE OPTIONAL ARGS, SO THIS MAY BE
% USEFUL.

% SEE BELOW FOR TEMPLATE FOR OBJECT CONSTRUCTORS

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

p_thr = .05;
dofdr = false;

dopartial = false;
dospearman = false;
dofigure = true;
doimage = false;
docircles = true;
dotext = true;
colorlimit = [-1 1]; % for correlations
max_radius = [];     % set by default to range(colorlimit)/4 or 0.5 if empty

partitions = [];
partitioncolors = {};
partitionlabels = {};

% names
var_names = [];

% text options
text_x_offset = .15;
text_y_offset = 0; % -.35; % negative is down
text_fsize = 16;
text_nonsig_color = [.3 .3 .3];
text_sig_color = [0 0 0];

[n, k] = size(X);

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

allowable_inputs = {'var_names' 'p_thr' 'dospearman' 'dopartial' 'dofdr' 'dofigure' 'doimage' 'docircles' 'dotext' 'colorlimit' 'text_x_offset' 'text_y_offset' 'text_fsize' 'text_nonsig_color' 'text_sig_color' 'partitions' 'partitioncolors' 'partitionlabels'};

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(max_radius)
    
    max_radius = range(colorlimit)/4;
    
end

if ~isempty(partitions)
    npartitions = length(unique(partitions));
    
    if isempty(partitioncolors)
        partitioncolors = scn_standard_colors(npartitions + 2);
        partitioncolors([1 3]) = [];  % remove red and blue (confusing)
    end
end


% Correlation stats
% --------------------------------------------------

if k > 30 && docircles
    disp('Warning: plotting correlations as circles will be slow with large number of variables');
end

if dopartial
    if dospearman
        [r, p] = partialcorr(double(X), 'type', 'Spearman');
    else
        [r, p] = partialcorr(double(X));
    end
    
elseif dospearman
    % Full pairwise, Spearman
    [r, p] = corr(double(X), 'Type', 'Spearman');
else
    % Full pairwise
    [r, p] = corr(double(X));
end

if dofdr
    % FDR-correct p-values across unique elements of the matrix
    
    trilp = tril(p, -1);
    wh = logical(tril(ones(size(p)), -1)); % Select off-diagonal values (lower triangle)
    trilp = double(trilp(wh));             % Vectorize and enforce double
    trilp(trilp < 10*eps) = 10*eps;        % Avoid exactly zero values
    p_thr = FDR(trilp, 0.05);
    
    if isempty(p_thr), p_thr = -Inf; end
    
end

is_sig = p < p_thr;

OUT.r = r;
OUT.p = p;
OUT.is_sig = is_sig;
OUT.p_thr = p_thr;
OUT.dospearman = dospearman;

% Names
% --------------------------------------------------
if ~isempty(var_names)
    var_names = strrep(var_names, '_', '');
end

OUT.var_names = var_names;

if dofigure
    create_figure('plotmatrix');
end

im_handle = imagesc(r, colorlimit);
set(gca, 'YDir', 'reverse', 'YTick', 1:k,  'YTickLabel', var_names, 'XTick', 1:k, 'XTickLabel', var_names, 'XTickLabelRotation', 45);


% Initial image matrix and colorbar
% --------------------------------------------------

colorbar
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm)

axis tight

% Get colormap for circles
r_ref_indx = linspace(-1, 1, size(cm, 1));

if ~doimage
    delete(im_handle)
end

% Partition colors
% --------------------------------------------------
if ~isempty(partitions)
    
    u = unique(partitions);
    if ~iscolumn(partitions), partitions = partitions'; end
    
    st = find(diff([0; partitions])) - 0.5;           % starting and ending values for each partition
    en = find(diff([partitions; npartitions + 1])) + 0.5;
    
    for i = 1:npartitions
        
        h1 = drawbox(st(i), en(i) - st(i), st(i), en(i) - st(i), 'k');
        set(h1, 'FaceColor', 'none', 'EdgeColor', 'k', 'lineWidth', 2);
        
        OUT.partition_handles(i) = h1;
        
        % partition color bars
        h1 = drawbox(st(i), en(i) - st(i), -0.75, 0.5, 'k');
        set(h1, 'EdgeColor', 'none', 'FaceColor', partitioncolors{i});
        OUT.partition_handles_hbars(i) = h1;
        
        h1 = drawbox(-0.75, 0.5, st(i), en(i) - st(i), 'k');
        set(h1, 'EdgeColor', 'none', 'FaceColor', partitioncolors{i});
        OUT.partition_handles_vbars(i) = h1;
        
        if ~isempty(partitionlabels)
            OUT.partition_handles_labels(i) = text(st(i), -0.55, partitionlabels{i}, 'FontSize', 14);
        end
        
    end
    
end


%% Draw circles
% --------------------------------------------------
% radius proportional to r, area proportional to r^2
% max radius = 0.5 to fit in unit square

if docircles
    
    for rr = 1:k
        
        for cc = 1:k
            
            % get color
            myr = r(rr, cc);
            [mymin, r_ind] = min((myr - r_ref_indx) .^ 2); % find closest
            mycolor = cm(r_ind, :);
            
            % draw filled circle
            [~, fhan] = circle([rr cc], max_radius * abs(myr), 'fill', mycolor);
            set(fhan, 'EdgeColor', 'none');
            
            % if significant, draw thicker border
            if rr == cc || is_sig(rr, cc)
                
                han = circle([rr cc], max_radius * abs(myr));
                set(han, 'LineWidth', 2, 'Color', mycolor ./ 1.5);
                
            end
            
        end
        
    end
    
end


%% Draw text
% --------------------------------------------------

if dotext
    
    for rr = 1:k
        
        for cc = 1:k
            
            % text
            myr = r(rr, cc);
            
            % if significant, bold weight
            if rr == cc || is_sig(rr, cc)
                
                text_han(rr, cc) = text(rr - text_x_offset, cc - text_y_offset, sprintf('%3.2f', myr), 'FontSize', text_fsize, 'FontWeight', 'b', 'Color', text_sig_color);
                
            else
                
                text_han(rr, cc) = text(rr - text_x_offset, cc - text_y_offset, sprintf('%3.2f', myr), 'FontSize', text_fsize, 'Color', text_nonsig_color);
                
            end
            
        end
        
    end
    
    OUT.text_han = text_han;
    
end

end % main function
