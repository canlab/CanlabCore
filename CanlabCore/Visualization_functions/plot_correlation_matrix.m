function OUT = plot_correlation_matrix(X, varargin)
% Plots an image of a correlation matrix with circles, text, image, or combinations
%
% :Usage:
% ::
%
%     OUT = plot_correlation_matrix(X, [optional inputs])
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
%       These include:
%       dofigure, doimage, docircles, dotext, colorlimit, text_x_offset, text_y_offset,  
%       text_fsize, text_nonsig_color, text_sig_color
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
% :References:
%   CITATION(s) HERE
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
dofigure = true;
doimage = false;
docircles = true;
dotext = true;
colorlimit = [-1 1]; % for correlations
max_radius = [];     % set by default to range(colorlimit)/4 or 0.5 if empty

% names
var_names = [];

% text options
text_x_offset = .15;
text_y_offset = 0; % -.35; % negative is down
text_fsize = 16;
text_nonsig_color = [.3 .3 .3];
text_sig_color = [0 0 0];

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

allowable_inputs = {'var_names' 'dofigure' 'doimage' 'docircles' 'dotext' 'colorlimit' 'text_x_offset' 'text_y_offset' 'text_fsize' 'text_nonsig_color' 'text_sig_color'};

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

% X = [D.LongShort D.PA_exp_pre D.PA_exp_post D.PA_beh_pre D.PA_beh_post];
% var_names = {'LongShort' 'PA_exp_pre' 'PA_exp_post' 'PA_beh_pre' 'PA_beh_post'};

% Correlation stats
% --------------------------------------------------

[r, p] = corr(X);
k = size(r, 1);

is_sig = p < p_thr;

OUT.r = r;
OUT.p = p;
OUT.is_sig = is_sig;
OUT.p_thr = p_thr;

% Names
% --------------------------------------------------
var_names = strrep(var_names, '_', '');

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
