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
%        Correlation matrix will be calculated, as per options specified
%        UNLESS:
%           X can be a structure compatible with output of ttest3d. 
%           It must have the fields 'r', 'p', and 'sig', each with a k x k matrix
%           r = correlation values
%           p = p-values
%           sig = signficance matrix (logical)
%        In this case, plot_correlation_matrix will ignore all calculation-related inputs
%
% :Optional Inputs:
%   **var_names:**
%        Followed by cell array of variable names
%
%   **p_thr:**
%        Followed by P-value threshold for significant correlations
%
%   **[OTHERS]:**
%
%   'spearman', 'Spearman', 'rank', 'dorank'   dospearman = true; 
%   'partial', 'Partial', 'partialcorr'        dopartial = true; 
%   'image'                                    doimage = true; docircles = false; 
%   'circles'                                  docircles = true; doimage = false;    
%   'notext'                                   dotext = false;
%   'fdr' 'FDR'                                dofdr = true; 
%   'nocalc' 'input_rmatrix'                   docalc = false; 
%                
%   'names', 'labels'   var_names = varargin{i+1}; varargin{i+1} = [];
%   Partitions: Note: partition labels MUST be sorted (ascending sequential) for this to work right.
%   'partitions'                                followed by k-length integer vector of partitions to plot with color bars
%   'partitioncolors'                           followed by cell array of strings for colors
%   'partitionlabels'                           followed by cell array of labels for each partition
%                
%       See code for other optional inputs controlling display
%       Enter a keyword followed by a value (e.g., true / false)
%       These include:
%       dofigure, dospearman, doimage, docircles, dotext, colorlimit, text_x_offset, text_y_offset,  
%       text_fsize, text_nonsig_color, text_sig_color
%       'dopartitions', followed by k-length integer vector defining partitions of variables  [optional] 
%       'partitioncolors', m-length cell vector defining rgb color triplets for each partition [optional] 
%       'p_thr', followed by p-value threshold
%
%       The default behavior is to plot circles with text values, UNLESS
%       your matrix is larger than 20 x 20, in which case it defaults to
%       image format with no labels.
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
%    Tor Wager: Aug 2019: Note: some options harmonize output with output structure of ttest3d
% ..

% BELOW IS A STANDARD TEMPLATE FOR DEFINING VARIABLE (OPTIONAL) INPUT
% ARGUMENTS. MANY FUNCTIONS NEED TO PARSE OPTIONAL ARGS, SO THIS MAY BE
% USEFUL.

% SEE BELOW FOR TEMPLATE FOR OBJECT CONSTRUCTORS

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

skip_calculation = false;  % not a valid input option. Used if structure with .r .sig .p entered instead of X

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
var_names = {};

% text options
text_x_offset = .15;
text_y_offset = 0; % -.35; % negative is down
text_fsize = 16;
text_nonsig_color = [.3 .3 .3];
text_sig_color = [0 0 0];

[n, k] = size(X);

if isstruct(X)  % Adjust if struct input
    [n, k] = size(X.r);
end

% adjust defaults if needed (these will be overridden by inputs below)
if k > 15
    dotext = false;
end

if k > 50
    docircles = false;
    doimage = true;
end

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% optional inputs with default values - each keyword entered will create a variable of the same name

allowable_inputs = {'var_names' 'p_thr' 'dospearman' 'dopartial' 'dofdr' 'dofigure' 'doimage' 'docircles' 'dotext' 'colorlimit' 'text_x_offset' 'text_y_offset' 'text_fsize' 'text_nonsig_color' 'text_sig_color' 'partitions' 'partitioncolors' 'partitionlabels'};

special_commands = {'spearman', 'Spearman', 'rank', 'dorank' 'partial', 'Partial', 'partialcorr' 'image' 'circles' 'notext' 'fdr' 'FDR' 'nocalc' 'input_rmatrix' 'nofigure' 'names', 'labels'};

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            case special_commands
                % do nothing; used later
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Allow for mapping of additional keywords. This harmonizes allowable
% inputs with other functions and makes entry of some options more
% intuitive.

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'spearman', 'Spearman', 'rank', 'dorank'}, dospearman = true; 
            case {'partial', 'Partial', 'partialcorr'},      dopartial = true; 
            case {'image'},                                  doimage = true; docircles = false; 
            case {'circles'},                                docircles = true; doimage = false;    
            case {'notext'},                                 dotext = false;
            case {'fdr' 'FDR'},                              dofdr = true; 
            case {'nocalc' 'input_rmatrix'},                 docalc = false; 
            case {'nofigure'},                               dofigure = false; hold on;
                
            case {'names', 'labels'}, var_names = varargin{i+1}; varargin{i+1} = [];
                
            case allowable_inputs   % do nothing
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Process variables that depend on values optional inputs: 

if isstruct(X)
    % X can be a structure compatible with output of ttest3d. 
    % It must have the fields 'r', 'p', and 'sig', each with a k x k matrix
    % r = correlation values
    % p = p-values
    % sig = signficance matrix (logical)
    
    r = X.r;
    p = X.p;
    sig = X.sig;
    k = size(r, 2);
    
    X = X.r;        % placeholder for attribute validation
    
    skip_calculation = true;
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

% -------------------------------------------------------------------------
% VALIDATE ATTRIBUTES OF INPUTS
% -------------------------------------------------------------------------

validateattributes(X,{'numeric'},{'2d'},'plot_correlation_matrix','X', 1);

logical_args = {'dofdr' 'false' 'dospearman' 'dofigure' 'doimage' 'docircles' 'dotext'};
for i = 1:length(logical_args)
    
    my_arg = eval([logical_args{i} ';']);
    validateattributes(my_arg,{'logical'},{'scalar'},'plot_correlation_matrix',logical_args{i});

end

cell_args = {'partitioncolors' 'partitionlabels' 'var_names'};
for i = 1:length(cell_args)
    
    my_arg = eval([cell_args{i} ';']);
    validateattributes(my_arg,{'cell'},{},'plot_correlation_matrix',cell_args{i});

end

vector_args = {'text_sig_color' 'text_nonsig_color'};
for i = 1:length(vector_args)
    
    my_arg = eval([vector_args{i} ';']);
    validateattributes(my_arg,{'numeric' 'vector' 'nonnegative' 'nonnan'},{},'plot_correlation_matrix',vector_args{i});

end



% Correlation stats
% --------------------------------------------------

if k > 50 && docircles
    disp('Warning: plotting correlations as circles will be slow with large number of variables');
end

% Remove nan values row-wise - 'rows' 'complete' method.
[wasnan, X] = nanremove(X);

if skip_calculation
% Do nothing - we already have r, p, sig

elseif dopartial
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

if ~skip_calculation && dofdr
    % FDR-correct p-values across unique elements of the matrix
    
    trilp = tril(p, -1);
    wh = logical(tril(ones(size(p)), -1)); % Select off-diagonal values (lower triangle)
    trilp = double(trilp(wh));             % Vectorize and enforce double
    trilp(trilp < 10*eps) = 10*eps;        % Avoid exactly zero values
    fdrthr = FDR(trilp, 0.05);
    
    if isempty(fdrthr), fdrthr = -Inf; end
    
    sig = p < fdrthr;
    
    % FDR-corrected - add output

    OUT.fdrsig = sig;
    OUT.fdrthr = fdrthr;
    
elseif ~skip_calculation
    % Uncorrected - add output
    
    sig = p < p_thr; 
    OUT.p_thr = p_thr;
    OUT.sigu = sig;

end

% Note: some options harmonize output with output structure of ttest3d
OUT.r = r;
OUT.p = p;
OUT.sig = sig;

if ~skip_calculation
    
    OUT.dospearman = dospearman;
    OUT.dopartial = dopartial;
    
end

% Names
% --------------------------------------------------
if ~isempty(var_names)
    var_names = strrep(var_names, '_', '');
end

OUT.var_names = var_names;

if dofigure
    create_figure('plotmatrix');
else
    hold on;
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
    
    boxwid = k/30;
    boxstart  = 0 - boxwid - .25;
    textstart = boxstart + .3 * boxwid;
    
    for i = 1:npartitions
        
        h1 = drawbox(st(i), en(i) - st(i), st(i), en(i) - st(i), 'k');
        set(h1, 'FaceColor', 'none', 'EdgeColor', 'k', 'lineWidth', 2);
        
        OUT.partition_handles(i) = h1;
        
        % partition color bars
        h1 = drawbox(st(i), en(i) - st(i), boxstart, boxwid, 'k');
        set(h1, 'EdgeColor', 'none', 'FaceColor', partitioncolors{i});
        OUT.partition_handles_hbars(i) = h1;
        
        h1 = drawbox(boxstart, boxwid, st(i), en(i) - st(i), 'k');
        set(h1, 'EdgeColor', 'none', 'FaceColor', partitioncolors{i});
        OUT.partition_handles_vbars(i) = h1;
        
        if ~isempty(partitionlabels)
            OUT.partition_handles_labels(i) = text(st(i), textstart, format_strings_for_legend(partitionlabels{i}), 'FontSize', 14);
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
            if rr == cc || sig(rr, cc)
                
                han = circle([rr cc], max_radius * abs(myr));
                set(han, 'LineWidth', 2, 'Color', mycolor ./ 1.5);
                
            end
            
        end
        
    end
    
end


%% Draw text
% --------------------------------------------------

text_han = [];

if dotext
    
    for rr = 1:k
        
        for cc = 1:k
            
            % text
            myr = r(rr, cc);
            
            % if significant, bold weight
            if rr == cc || sig(rr, cc)
                
                text_han(rr, cc) = text(rr - text_x_offset, cc - text_y_offset, sprintf('%3.2f', myr), 'FontSize', text_fsize, 'FontWeight', 'b', 'Color', text_sig_color);
                
            else
                
                text_han(rr, cc) = text(rr - text_x_offset, cc - text_y_offset, sprintf('%3.2f', myr), 'FontSize', text_fsize, 'Color', text_nonsig_color);
                
            end
            
        end
        
    end
    
    OUT.text_han = text_han;
    
end

end % main function
