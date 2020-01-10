function [han, X, Y, slope_stats] = line_plot_multisubject(X, Y, varargin)
% Plots a scatterplot with multi-subject data, with one line per subject
% in a unique color.
%
% :Usage:
% ::
%
%    [han, X, Y, slope_stats] = line_plot_multisubject(X, Y, varargin)
%
% :Inputs:
%
%   **X and Y:**
%        are cell arrays, one cell per upper level unit (subject)
%
% varargin:
%
%   **'n_bins:**
%        pass in the number of point "bins".  Will divide each subj's trials
%        into bins, get the avg X and Y per bin, and plot those points.  
%
%   **'noind':**
%        suppress points
%
%   **'subjid':**
%        followed by integer vector of subject ID numbers. Use when
%        passing in vectors (with subjects concatenated) rather than
%        cell arrays in X and Y
%
%   **'center':**
%        subtract means of each subject before plotting
%
%   **'colors':**
%        followed by cell array length N of desired colors, rgb specification,
%        for each line.  if not passed in, will use custom_colors.
%        Group average lines use color{1} for line and color{2} for fill.
%
%   **'gcolors':***
%        Group average line colors, {[r g b] [r g b]} for line and point
%        fill, respectively
%
%   **'MarkerTypes':**
%        followed by char string.  if not passed in, uses
%        'osvd^<>ph' by default
%
%   **'group_avg_ref_line':**
%        will make a reference line for the group avg
%
% :Outputs:
%
%   **han:**
%        handles to points and lines
%
%   **X, Y:**
%        new variables (binned if bins requested)
%
%   **slope_stats** 
%        slope of linear relationship for each person
%         
% :Examples:
% ::
%
%    % Complete example with data generation and results
%    % -----------------------------------------------------------
%    % Create data for 5 simulated subjects, 10 observations each, random intercept, random positive slope:
%    for i = 1:5, X{i} = randn(10, 1); Y{i} = rand(1) * X{i} + .3 * randn(10, 1) + randn(1); end
%    han = line_plot_multisubject(X, Y)
% %
%    % Plot the results three ways, with different data transformations:
%     create_figure('Line plot multisubject', 1, 3);
%     han = line_plot_multisubject(X, Y);
%     title('Raw data (no transformation');
%     subplot(1, 3, 2);
%     han = line_plot_multisubject(X, Y, 'center');
%     title('Centered within-person');
%     subplot(1, 3, 3);
%     han = line_plot_multisubject(X, Y, 'zscore');
%     title('Z-scored within-person');
%
%
%   % -----------------------------------------------------------
  % Example creating bins of data within-person, useful for many within-person observations

  % Create data for 5 simulated subjects, 100 observations each, random intercept, random positive slope:
%    for i = 1:5, expect{i} = randn(100, 1); pain{i} = rand(1) * expect{i} + .3 * randn(100, 1) + randn(1); end

%   % Plot with bins, custom colors and points:
%   create_figure('Line plot multisubject with bins');
   %[han, Xbin, Ybin] = line_plot_multisubject(expect, pain, 'n_bins', 4, 'group_avg_ref_line', 'MarkerTypes', 'o', 'colors', custom_colors([1 .7 .4], [1 .7 .4], 100));

%   % -----------------------------------------------------------
%
% Center within subjects and bin, then calculate correlation of
% within-subject variables:
% ::
%
%    create_figure('lines'); [han, Xbin, Ybin] = line_plot_multisubject(stats.Y, stats.yfit, 'n_bins', 7, 'center');
%    corr(cat(1, Xbin{:}), cat(1, Ybin{:}))

%    Programmer's notes:
%    12/22/19 - Marta: added z-scoring option, added and clarified
%    different flavors of r (overall r, within-subject r, between-subject r)
%    01/05/20 - Marta: added NaN handling for r_within for case when all input
%    values of a subject are 0)


% -------------------------------------------------------------------------
% Defaults and inputs
% -------------------------------------------------------------------------

docenter = false;
dozscore = false; 
doind = true;
dolines = true;
group_avg_ref_line = false;

for i=1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'n_bins'
                n_bins = varargin{i+1};
                %bin_size = length(X{1}) / n_bins;
                %if rem(length(X{1}), n_bins) ~= 0, error('Num trials must be divisible by num bins'), end;
                
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
                
            case {'center', 'docenter'}
                docenter = 1;
                
            case {'zscore', 'dozscore'}
                dozscore = 1;
                
            case 'colors'
                colors = varargin{i+1};
                
            case 'MarkerTypes'
                mtypes = varargin{i+1};
                
            case 'noind'
                doind=0;
                
            case 'nolines'
                dolines = 0;
                
            case 'group_avg_ref_line'
                group_avg_ref_line = 1;
                
            case 'gcolors' % Mj added!!
                gcolors = varargin{i+1};
        end
    end
end

if ~iscell(X) || ~iscell(Y)
    error('X and Y should be cell arrays, one cell per line to plot.');
end

N = length(X);

% Set color options and check
% -------------------------------------------------------------------------

if ~exist('colors', 'var') %,colors = scn_standard_colors(N); end
    colors = custom_colors([1 .5 .4], [.8 .8 .4], N);
end
if ~iscell(colors), error('Colors should be cell array of color specifications.'); end
if length(colors) < N, colors = repmat(colors, 1, N); end

if ~exist('mtypes', 'var'), mtypes = 'osvd^<>ph'; end

if ~exist('gcolors', 'var') 
    % MJ added.  Line then fill.
    gcolors = colors(1:2);
end
if ~iscell(gcolors), gcolors = {gcolors}; end
    
% replicate for fill color if 2nd color is empty
if length(gcolors) < 2, gcolors{2} = gcolors{1}; end

hold on

% -------------------------------------------------------------------------
% Center/scale data as requested
% -------------------------------------------------------------------------

% Marta 12/22/19 (+Tor):
% Calc mean and std per subject before doing anything else to data (excluding NaNs)
% This is useful for reporting on input data state
X_between=cellfun(@nanmean, X);
Y_between=cellfun(@nanmean, Y);

X_std=cellfun(@nanstd, X); % 01/05/2020 Marta corrected nanmean to nanstd
Y_std=cellfun(@nanstd, Y); 

% Save some text for a report to be printed later

X_is_centered = all(abs(X_between)) < 1e-10; % fix this to be valid for 0s too
Y_is_centered = all(abs(Y_between)) < 1e-10;
X_is_zscored = X_is_centered && all(abs(X_std - 1)) < 1e-10;
Y_is_zscored = Y_is_centered && all(abs(Y_std - 1)) < 1e-10;

if X_is_zscored, X_data_str = 'Z-scored'; elseif X_is_centered, X_data_str = 'Centered'; else X_data_str = 'No centering or z-scoring'; end
if Y_is_zscored, Y_data_str = 'Z-scored'; elseif Y_is_centered, Y_data_str = 'Centered'; else Y_data_str = 'No centering or z-scoring'; end

report_str{1} = 'Input data:';
report_str{2} = sprintf('X scaling: %s', X_data_str);
report_str{3} = sprintf('Y scaling: %s', Y_data_str);
report_str{4} = '\nTransformations:';
 
% zscore, if asked for (removes mean AND divides by std) 
if dozscore
    X = cellfun(@scale, X, 'UniformOutput', false);
    Y = cellfun(@scale, Y, 'UniformOutput', false);
   report_str{5} = 'X and Y Z-scored before plot';

elseif docenter
% center, if asked for (removes mean)
    X = cellfun(@(x) scale(x, 1), X,  'UniformOutput', false);
    Y = cellfun(@(x) scale(x, 1), Y,  'UniformOutput', false);
    report_str{5} = 'X and Y centered (forced mean-zero) before plot';

else
   report_str{5} = 'No data transformations before plot';
end


% -------------------------------------------------------------------------
% Plot points and lines
% -------------------------------------------------------------------------

for i = 1:N
    
    % choose marker
    whm = min(i, mod(i, length(mtypes)) + 1);
    
    if length(X{i}) == 0 || all(isnan(X{i})) ||  all(isnan(Y{i}))
        % empty
        continue
        
    elseif length(X{i}) ~= length(Y{i})
        error(['Subject ' num2str(i) ' has unequal elements in X{i} and Y{i}. Check data.'])
    end

    % plot points in bins
    if exist('n_bins', 'var')
        if n_bins ~= 0
            points = zeros(n_bins,2);
            t = sortrows([X{i} Y{i}],1); % not really even needed
            
            x = t(:, 1);
            bins = prctile(x, linspace(0, 100, n_bins + 1));
            bins(end) = Inf;
            
            for j=1:n_bins %make the bins
                
                wh = x >= bins(j) & x < bins(j+1);
                
                points(j,:) = nanmean(t(wh, :));
                
                ste_points(j, :) = nanstd(t(wh, :)) ./ sqrt(sum(wh));
                
                %points(j,:) = nanmean(t( bin_size*(j-1)+1 : bin_size*j ,:));
            end
            
            X{i} = points(:, 1);
            Y{i} = points(:, 2);
            
            %han.point_handles(i) = plot(points(:,1), points(:,2), ['k' mtypes(whm(1))], 'MarkerFaceColor', colors{i}, 'Color', max([0 0 0; colors{i}.*.7]));
        end
    end
    
    % plot ref line
    b(i,:) = glmfit(X{i}, Y{i});
    
    if dolines
        han.line_handles(i) = plot([min(X{i}) max(X{i})], [b(i,1)+b(i,2)*min(X{i}) b(i,1)+b(i,2)*max(X{i})], 'Color', colors{i}, 'LineWidth', 1);
    end
    
    % plot all the points
    if doind
        han.point_handles(i) = plot(X{i}, Y{i}, ['k' mtypes(whm(1))], 'MarkerSize', 3, 'MarkerFaceColor', colors{i}, 'Color', max([0 0 0; colors{i}]));
    end
end % subject loop

% Stats on slope
% -----------------------------------------------------------------
slope_stats.b = b;
[~, slope_stats.p, ~, tmpstat] = ttest(b(:, 2));
slope_stats.t = tmpstat.tstat;
slope_stats.df = tmpstat.df;
    
Xc = cat(1, X{:});
Yc = cat(1, Y{:});
[wasnan, Xc, Yc] = nanremove(Xc, Yc);
slope_stats.r = corr(Xc, Yc);
slope_stats.wasnan = wasnan;

% Marta 12/22/19 

% Overall r 
% -----------------------------------------------------------------
% Clarify consequences of removing mean and z-scoring
if docenter
    report_str{6} = sprintf('\nCorrelations:\nr = %3.2f across all observations, removing subject mean (X and Y are centered)', slope_stats.r);
elseif dozscore
    report_str{6} = sprintf('\nCorrelations:\nr = %3.2f across all observations, based on Z-scored X and Y data', slope_stats.r);
else
    report_str{6} = sprintf('\nCorrelations:\nr = %3.2f across all observations, based on untransformed input data', slope_stats.r);
end


% Within-person r (each subject's X and Y correlated, then the r values averaged across subjects) 
% -----------------------------------------------------------------
[wasnan, Xr]= cellfun(@nanremove, X, 'UniformOutput', false); % Marta 1/7/20 unique name for nan-removed cell arrays 
[wasnan, Yr]= cellfun(@nanremove, Y, 'UniformOutput', false);
r_within = cellfun(@corr, Xr,Yr);

% Marta 01/05/2020
if (nnz(isnan(r_within)))>0;
    sprintf('%3.0f subject(s) have an r value of NaN, will remove remove NaNs', ...
       nnz(isnan(r_within)))
    r_within=r_within(~isnan(r_within));
end
    
slope_stats.r_within = mean(r_within);
slope_stats.r_within_std = std(r_within);

report_str{8} = sprintf('Average within-person r = %3.2f +- %3.2f (std)', ...
    slope_stats.r_within, slope_stats.r_within_std);

if dozscore
    report_str{9} = sprintf('* Note that the overall r and average within-person r are the same because X and Y data are z-scored\n');

elseif docenter
    report_str{9} = sprintf('* Note that the overall r and average within-person r may be similar because subject mean is removed\n');  

else
    report_str{9} = '';
end 

% Between-person r
% -----------------------------------------------------------------
slope_stats.r_between = corr(X_between', Y_between');
report_str{10} = sprintf('Between-person r (across subject means) = %3.2f', ...
    slope_stats.r_between);


% Stats on slopes across subjects
% -----------------------------------------------------------------
report_str{7} = sprintf('\nStats on slopes after transormation, subject is random effect: \nMean b = %3.2f, t(%3.0f) = %3.2f, p = %3.6f, num. missing: %3.0f\n', ...
    nanmean(slope_stats.b(:, 2)), slope_stats.df, slope_stats.t, slope_stats.p, sum(slope_stats.wasnan));


% Print report
% -----------------------------------------------------------------
% All stats/text have been collected now

canlab_print_legend_text(report_str{:});


% Individual points
% -----------------------------------------------------------------
% get rid of bad handles for missing subjects
if doind
    han.point_handles = han.point_handles(ishandle(han.point_handles) & han.point_handles ~= 0);
end

if dolines
    han.line_handles = han.line_handles(ishandle(han.line_handles) & han.line_handles ~= 0);
end

%the correlation
%r=corr(cat(1,detrend(X{:},'constant')), cat(1,detrend(Y{:}, 'constant')), 'rows', 'complete')

% plot the average ref line
% -----------------------------------------------------------------
if group_avg_ref_line
    
    if exist('n_bins', 'var') && n_bins ~= 0
        % PLOT GROUP BIN WITH CROSSHAIR STD. ERRORS FOR EACH BIN
        
        XX = cat(2, X{:})';
        YY = cat(2, Y{:})';
        
        % means and standard errors
        mX = nanmean(XX);
        sX = ste(XX);
        mY = nanmean(YY);
        sY = ste(YY);
        
        h = sepplot(mX, mY, .7, 'color', gcolors{1}, 'linewidth', 4);
        
        h2 = errorbar_horizontal(mX, mY, sX, 'o', 'color', gcolors{1}, 'linewidth', 3, 'markersize', 8, 'markerfacecolor', gcolors{2});
        h1 = errorbar(mX, mY, sY, 'o', 'color', gcolors{1}, 'linewidth', 3, 'markersize', 8, 'markerfacecolor', gcolors{2});
        set(h2, 'linewidth', 3, 'color', gcolors{1});
        
        han.grpline_handle = h;
        han.grpline_err_handle = [h1 h2];
        else
            
    avg_b = mean(b);
    Xs = cat(2,X{:});
    minX = prctile(Xs(1,:),5);
    maxX = prctile(Xs(end,:),95);
    han.grpline_handle = plot([minX maxX], [avg_b(1)+avg_b(2)*minX avg_b(1)+avg_b(2)*maxX], 'Color', 'k', 'LineWidth', 6);
end

end % function
