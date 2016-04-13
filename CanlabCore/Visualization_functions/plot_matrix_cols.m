function han = plot_matrix_cols(X, varargin)
% :Usage:
% ::
%
%    han = plot_matrix_cols(X, [method], [x-values], [colors cell], [linewidth], [axis limits])
%
% Plot line plots showing each column of a matrix as a vertical or
% horizontal line
% ::
%
%    han = plot_matrix_cols(X)
%
% :Optional Inputs:
%
%   **method:**
%        plot_matrix_cols(X, 'horiz')
%
%        plot_matrix_cols(X, 'vertical')
%
%   **x-values:**
%        plot_matrix_cols(X, 'horiz', x_in_secs)
%
%   **colors:**
%        plot_matrix_cols(X, 'horiz', [], {'r' 'g' 'b' 'y' 'm'})
%
%   **Linewidth:**
%        plot_matrix_cols(X, 'horiz', [], {'r' 'g' 'b' 'y' 'm'}, 2)
%
%   **axis limits:**
%        plot_matrix_cols(X, 'horiz', [], [], [], [0 10])
% 
% legacy 'method' string:
% within denoising, plot_matrix_cols(X, 'denoising') to make red plots
%
% ..
%    tor wager, feb 07
%    updated 8/2015 by tor to add colors, make faster
% ..

meth = 'horiz'; % Default inputs
lineWidth = 1;

% Calculated default variables
% -----------------------------------------------------------

[t, k] = size(X);
colors = repmat({[.5 .5 .5]}, 1, k);
tvec = 1:t;
axislimits = [0 k+1];

% Optional inputs
% -----------------------------------------------------------
if ~isempty(varargin)
    
    if ~isempty(varargin{1}), meth = varargin{1}; end
    
    if length(varargin) > 1 && ~isempty(varargin{2})
        tvec = varargin{2};
    end
    
    if length(varargin) > 2 && ~isempty(varargin{3})
        colors = varargin{3};
    end
    
    if length(varargin) > 3 && ~isempty(varargin{4})
        lineWidth = varargin{4};
    end
    
    if length(varargin) > 4 && ~isempty(varargin{5})
        axislimits = varargin{5};
    end
    
end

% scale vals in each col of X to range = ~0.9
% -----------------------------------------------------------
X = scale(X, 1);
rng = range(X);
X = X ./ (1.2 * rng(ones(t, 1), :));


% Set axes
% -----------------------------------------------------------
hold on;

switch meth
    
    case {'horiz', 'horizontal', 'h'}
        set(gca, 'YLim', [0 k+1]);
        
    case {'vert', 'vertical', 'v'}
        set(gca, 'XLim', [0 k+1]);
        
end

set(gca, 'YDir', 'Reverse');

% Loop and plot
% -----------------------------------------------------------
for i = 1:k
    
    switch meth
        case {'horiz', 'horizontal', 'h'}
            
            %plot(tvec, i, 'Color', [.7 .7 .7]);
            % negative b/c y dir is reversed
            han(i) = plot(tvec, -X(:, i) + i, '-', 'Color', colors{i}, 'LineWidth', lineWidth);
            
            
            
        case {'vert', 'vertical', 'v'}
            
            %plot(i, tvec, 'Color', [.7 .7 .7]);
            
            han(i) = plot(X(:, i) + i, tvec, 'Color', colors{i}, 'LineWidth', lineWidth);
            
            
        case {'denoising'}
            plot(tvec, i, 'Color', [1 0 0]);
            % negative b/c y dir is reversed
            han(i) = plot(tvec, -X(:, i) + i, 'r', 'LineWidth', lineWidth);
            
        otherwise
            error('Unknown method.');
    end
    
    %drawnow
    
end


end
