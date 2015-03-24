function [han, X, Y] = line_plot_multisubject(X, Y, varargin)
% Plots a scatterplot with multi-subject data, with one line per subject
% in a unique color.
%
% [han, X, Y] = line_plot_multisubject(X, Y, varargin)
%
% Inputs:
% % X and Y are cell arrays, one cell per upper level unit (subject)
%
% varargin:
%    'n_bins' = pass in the number of point "bins".  Will divide each subj's trials
%               into bins, get the avg X and Y per bin, and plot those
%               points.  
%    'noind' - suppress points
%
%    'subjid' - followed by integer vector of subject ID numbers. Use when
%    passing in vectors (with subjects concatenated) rather than cell arrays in X and Y
%
%    'center' - subtract means of each subject before plotting
%    'colors' - followed by array size N of desired colors.  if not passed
%    in, will use scn_standard_colors
%    'MarkerTypes' - followed by char string.  if not passed in, uses
%    'osvd^<>ph' by default
%
%    'group_avg_ref_line' - will make a reference line for the group avg
%
% Outputs:
% han = handles to points and lines
% X, Y = new variables (binned if bins requested)
%
% Examples:
% -------------------------------------------------------------------------
% for i = 1:20, X{i} = randn(4, 1); Y{i} = X{i} + .3*randn(4, 1) + randn(1); end
% han = line_plot_multisubject(X, Y)
%
% Center within subjects and bin, then calculate correlation of
% within-subject variables:
% create_figure('lines'); [han, Xbin, Ybin] = line_plot_multisubject(stats.Y, stats.yfit, 'n_bins', 7, 'center');
% corr(cat(1, Xbin{:}), cat(1, Ybin{:}))

docenter = 0;
doind = 1;
group_avg_ref_line = 0;

for i=1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'n_bins'
                n_bins = varargin{i+1};
                bin_size = length(X{1}) / n_bins;
                if rem(length(X{1}), n_bins) ~= 0, error('Num trials must be divisible by num bins'), end;
                
            case 'subjid'
                subjid = varargin{i + 1};
                if iscell(X) || iscell(Y)
                    error('X and Y should be vectors when using subjid, not cell arrays.');
                end
                u = unique(subjid);
                for i = 1:length(u)
                    XX{i} = X(subjid == u(i));
                    YY{i} = Y(subjid == u(i));
                end
                X = XX;
                Y = YY;
                
            case 'center'
                docenter = 1;
            case 'colors'
                colors = varargin{i+1};
            case 'MarkerTypes'
                mtypes = varargin{i+1};
            case 'noind'
                doind=0;
            case 'group_avg_ref_line'
                group_avg_ref_line = 1;
        end
    end
end

if ~iscell(X) || ~iscell(Y)
    error('X and Y should be cell arrays, one cell per line to plot.');
end


N = length(X);

if ~exist('colors', 'var'), colors = scn_standard_colors(N); end
if ~exist('mtypes', 'var'), mtypes = 'osvd^<>ph'; end

hold on

for i = 1:N
    
    % choose marker
    whm = min(i, mod(i, length(mtypes)) + 1);
    
    if length(X{i}) == 0 || all(isnan(X{i})) ||  all(isnan(Y{i}))
        % empty
        continue
        
    elseif length(X{i}) ~= length(Y{i})
        error(['Subject ' num2str(i) ' has unequal elements in X{i} and Y{i}. Check data.'])
    end
    
    % centere, if asked for
    if docenter
        X{i} = scale(X{i}, 1);
        Y{i} = scale(Y{i}, 1);
    end
    
    
    % plot points in bins
    if exist('n_bins', 'var')
        if n_bins ~= 0
            points = zeros(n_bins,2);
            t = sortrows([X{i} Y{i}],1);
            for j=1:n_bins %make the bins
                points(j,:) = nanmean(t( bin_size*(j-1)+1 : bin_size*j ,:));
            end
            
            X{i} = points(:, 1);
            Y{i} = points(:, 2);
            
            %han.point_handles(i) = plot(points(:,1), points(:,2), ['k' mtypes(whm(1))], 'MarkerFaceColor', colors{i}, 'Color', max([0 0 0; colors{i}.*.7]));
        end
    end
    
    % plot ref line
    b(i,:) = glmfit(X{i}, Y{i});
    
    han.line_handles(i) = plot([min(X{i}) max(X{i})], [b(i,1)+b(i,2)*min(X{i}) b(i,1)+b(i,2)*max(X{i})], 'Color', colors{i}, 'LineWidth', 1);
    
    % plot all the points
    if doind
        han.point_handles(i) = plot(X{i}, Y{i}, ['k' mtypes(whm(1))], 'MarkerFaceColor', colors{i}, 'Color', max([0 0 0; colors{i}]));
    end
end % subject loop

% get rid of bad handles for missing subjects
if doind
    han.point_handles = han.point_handles(ishandle(han.point_handles) & han.point_handles ~= 0);
end
han.line_handles = han.line_handles(ishandle(han.line_handles) & han.line_handles ~= 0);

%the correlation
%r=corr(cat(1,detrend(X{:},'constant')), cat(1,detrend(Y{:}, 'constant')), 'rows', 'complete')

% plot the average ref line
if group_avg_ref_line
    avg_b = mean(b);
    Xs = cat(2,X{:});
    minX = prctile(Xs(1,:),5);
    maxX = prctile(Xs(end,:),95);
    plot([minX maxX], [avg_b(1)+avg_b(2)*minX avg_b(1)+avg_b(2)*maxX], 'Color', 'k', 'LineWidth', 6);
end

end % function
