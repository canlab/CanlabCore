function han = plot_matrix_cols(X, varargin)
%
% Plot line plots showing each column of a matrix as a vertical or
% horizontal line
%
% han = plot_matrix_cols(X)
% plot_matrix_cols(X, 'horiz')
% plot_matrix_cols(X, 'vertical')

% within denoising, plot_matrix_cols(X, 'denoising') to make red plots
%
% tor wager, feb 07

meth = 'horiz';
if ~isempty(varargin), meth = varargin{1}; end

[t, k] = size(X);

% scale vals in each col of X to range = 1
X = scale(X, 1);
rng = range(X);
X = X ./ rng(ones(t, 1), :);

% plot
hold on;


for i = 1:k

    switch meth
        case {'horiz', 'horizontal', 'h'}
            set(gca, 'YLim', [0 k+1]);
            plot(1:t, i, 'Color', [.5 .5 .5]);
            % negative b/c y dir is reversed
            han(i) = plot(1:t, -X(:, i) + i, 'k');

            

        case {'vert', 'vertical', 'v'}
            set(gca, 'XLim', [0 k+1]);
            plot(i, 1:t, 'Color', [.5 .5 .5]);
            
            han(i) = plot(X(:, i) + i, 1:t, 'k');

            
        case {'denoising'}
             plot(1:t, i, 'Color', [1 0 0]);
            % negative b/c y dir is reversed
            han(i) = plot(1:t, -X(:, i) + i, 'r');
            
        otherwise
            error('Unknown method.');
    end

    drawnow
    
end

set(gca, 'YDir', 'Reverse');

end