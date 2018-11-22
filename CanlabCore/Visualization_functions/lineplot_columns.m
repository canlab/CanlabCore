function out = lineplot_columns(dat_matrix, varargin)
% :Usage:
% ::
%
%    out = lineplot_columns(dat_matrix, varargin)
%
% :Default values: (see below for how to change)
%   - w = 3;                          % width
%   - color = 'k';                    % color
%   - wh = true(size(dat_matrix));    % which observations
%   - x = 1:size(dat_matrix, 2);      % x values
%   - marker = 'o';                   % marker style
%   - linestyle = '-';                % line
%   - markersize = 8;                 % markersize
%   - markerfacecolor = [.5 .5 .5];   % face color
%   - dowithinste = 0;                % enter 'within' to get within-ss ste
%   - atleast = 1;                    % 'atleast' followed by n for at least n valid obs to plot 
%   - doshading = 0;                  % shaded vs. error-bar plots
%   - CIs95 = 0;                      % 'CIs', '95CIs', or 'CIs95' for error bars to show 95% CIs. Not supported for within subj error
%
% :Inputs:
%
%   **dat_matrix:**
%        is usually a rectangular matrix with rows = observations,
%        columns = variables.
%
%        It can also be a cell array with column vectors in each cell, for unequal
%        numbers of observations, but then the rows will not be the same
%        observations across variables.
%
% :Optional Inputs: followed by values:
%
%   {'w', 'color', 'x', 'marker', 'linestyle', 'markersize', 'markerfacecolor', 'wh'}
%
% :Keywords:
%
% {'within', 'dowithinste'}, 'shade', 'atleast', 'CIs' <- for 95% CI's
%
% :Examples:
% ::
%
%    out = lineplot_columns(dat_matrix, 'w', 3, 'color', 'r', 'markerfacecolor', [1 .5 0], 'wh', ispain);
%    out = lineplot_columns(dat_matrix, 'w', 3, 'color', [0 1 0], 'markerfacecolor', [1 .5 0], 'within');
%    out = lineplot_columns(dat_matrix, 'w', 3, 'color', 'b', 'markerfacecolor', [0 .5 1], 'within');
%
%    % shaded error regions
%    out = lineplot_columns(hotopen, 'color', 'r', 'marker', 'none', 'w', 1, 'shade');

w = 3;                          % width
color = 'k';                    % color
wh = true(size(dat_matrix));    % which observations
x = 1:size(dat_matrix, 2);      % x values
marker = 'o';                   % marker style
linestyle = '-';                % line
markersize = 8;                 % markersize
markerfacecolor = [.5 .5 .5];   % face color
dowithinste = 0;                % enter 'within' to get within-ss ste
atleast = 1;
doshading = 0;                  % shaded vs. error-bar plots
CIs95 = 0;                      % 95% CI's.  Does not work with within subj error

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            % functional commands
            case {'w', 'color', 'x', 'marker', 'linestyle', 'markersize', 'markerfacecolor', 'wh'}
                str = [varargin{i} ' = varargin{i + 1};'];
                eval(str)
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case {'within', 'dowithinste'}
                dowithinste = 1;
                
            case 'atleast'
                atleast = varargin{i + 1};
                
            case 'shade'
                doshading = 1;
                
            case {'cis', '95cis', 'cis95'}
                CIs95 = 1;
                
            otherwise
                error('Unknown input.');
        end
    end
end

% If cells, pad with NaNs and make matrix
if iscell(dat_matrix)
    fprintf('Converting from cell array, padding with NaNs if necessary.\n');
    
    Csize = cellfun(@size, dat_matrix, 'UniformOutput', false);
    Csize = cat(1, Csize{:});
    [mx, wh] = max(Csize(:, 1));
    
    for i = 1:length(dat_matrix)
        dat_matrix{i} = padwithnan(dat_matrix{i}, dat_matrix{wh}, 1);
    end
    
    dat_matrix = cat(2, dat_matrix{:});
    
end

%%

dat_matrix(~wh) = NaN;

wh_omit = sum(~isnan(dat_matrix)) < atleast;

dat_matrix(:, wh_omit) = [];
x(wh_omit) = [];

if dowithinste
    
    [out.ste, stats] = barplot_get_within_ste(dat_matrix);
    out.ste = repmat(out.ste, 1, size(dat_matrix, 2));
    
else
    out.ste = ste(dat_matrix);
    
    ts = tinv(.975, size(dat_matrix, 1));
    out.CI95 = ts * out.ste;
end

out.m = nanmean(dat_matrix);

hold on;

if doshading
    if CIs95
        out.err_han = fill_around_line(out.m, out.CI95, color, x);
    else
        out.err_han = fill_around_line(out.m, out.ste, color, x);
    end
else
    if CIs95
        out.err_han = errorbar(x, out.m, out.CI95, 'Color', color, 'LineWidth', w);
    else    
        out.err_han = errorbar(x, out.m, out.ste, 'Color', color, 'LineWidth', w);
    end
    
    set(out.err_han, 'LineStyle', 'none')
end

out.line_han = plot(x, out.m, 'Color', color, 'Marker', marker, 'LineStyle', linestyle, ...
    'MarkerSize', markersize, 'MarkerFaceColor', markerfacecolor, 'LineWidth', w);

end % function

