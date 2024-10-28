function [est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(dat, varargin)
% Find outliers based on a combination of rmssd (DVARS), robust spatial
% variance (MAD), Mahalanobis distance based on covariance and correlation
% matrices, and images with >25% missing values
%
% :Usage:
% ::
%
% [est_outliers_uncorr, est_outliers_corr, outlier_table] = outliers(dat, ['noverbose', 'notimeseries'])
%
% Images usually change slowly over time, and sudden changes in intensity can also often be a sign of bad things
% -- head movement artifact or gradient misfires, interacting with the magnetic field to create distortion
% across the brain.
%
% RMSSD (DVARS) tracks large changes across successive images, regardless of what the sign of the changes is or where they are.
% In addition, images with unusually high spatial standard deviation across voxels may be outliers with image
% intensity distortions in some areas of the image but not others (e.g., bottom half of brain vs. top half,
% or odd vs. even slices).
% This matrix can be added to your design matrix as a set of nuisance covariates of no interest.
%
% :Optional Inputs:
%
%   **'madlim'**
%   Limit for median absolute deviation (MAD) for rmssd and spatial abs. deviation
%
%   **'noverbose'**
%   Suppress verbose output
%
%   **'noplot'**
%   Suppress plot output
%
%   **'fullplot'**
%   A more detailed plot of each criterion
%
%   **'notimeseries'**
%   Suppress time series-specific metrics -- use for 2nd-level contrasts or beta series
%
%   **'fd'**
%   called framewise_displacement() to generate framewise_displacement
%   indicators. Requires movement matrix to be passed in as an argument.
%
%
% :Outputs:
%
%   **est_outliers_uncorr**
%   Logical vector of outliers at uncorrected thresholds
%
%   **est_outliers_corr**
%   Logical vector of outliers at corrected thresholds (more conservative)
%
%   **outlier_tables**
%   Tables of outliers under various criteria, and indicator matrices for
%   nuisance regressors
%
%     outlier_indicator_table:          Table of logical indicator vectors for each criterion measure (Matlab Table object)
%     score_table:                      Table of criterion scores (Matlab Table object)
%     summary_table:                    Table of outlier counts for each measure, and overall
%     est_outliers_uncorr               Logical vector of outliers at uncorrected thresholds
%     est_outliers_corr                 Logical vector of outliers at corrected thresholds
%     outlier_regressor_matrix_uncorr:  Indicator matrix for nuisance regressors, uncorrected
%     outlier_regressor_matrix_corr:    Indicator matrix for nuisance regressors, uncorrected
%
% :Examples:
% ::
% % -------------------------------------------------------------------------
% % Load a multi-study dataset, rescale it, and identify/plot outliers
% % Use 'notimeseries' option because this is not a time series dataset
% 
% obj = load_image_set('kragel18_alldata');
% obj2 = rescale(obj, 'l2norm_images');     % normalize heterogeneous datasets
% [est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(obj2, 'notimeseries');
%
% % -------------------------------------------------------------------------
%
% % -------------------------------------------------------------------------
% % Load a dataset from CANlab 2nd-level batch script output and assess outliers
% % across condition images. Select subjects with no outliers in any
% % condition
% 
% [est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(obj2, 'notimeseries');
% load('data_objects.mat')  % Load DATA_OBJ
% load('image_names_and_setup.mat'); % Load DAT
% obj = cat(DATA_OBJ{:});
% [est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(obj, 'notimeseries');
% 
% % -------------------------------------------------------------------------

% Programmers Notes:
% ..
%    Created 11/6/2021 by Tor Wager, from a combination of other code
%    (default: based on estimated outliers at 3 standard deviations.)
%
%    Updated 10/1/2024 by Michael Sun PhD to add fd argument for
%    framewise_displacement.
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------
madlim = 3;             % Also Z-score limit for global means
dotimeseries = true;    % Adds rmssd/dvars, time series-specific outliers
doverbose = true;
verbosestr = 'doverbose';
doplot = true;
dobriefplot = true;
dofd = false;

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'madlim' 'doverbose' 'dotimeseries' 'plot', 'fd'};

keyword_inputs = {'noverbose' 'notimeseries' 'noplot' 'fullplot'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                if strcmp(varargin{i}, 'fd')
                    dofd = true;
                end
                
            case keyword_inputs
                % Skip, deal with these below
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% 2nd pass: Keyword inputs. These supersede earlier inputs
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'noverbose'}
                doverbose = false;
                
            case {'notimeseries'}
                dotimeseries = false;
                
            case {'noplot'}
                doplot = false;

            case 'fullplot'
                dobriefplot = false;

        end
    end
end

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Global mean and variance-related
% -------------------------------------------------------------------------

if doverbose
    disp('______________________________________________________________')
    disp('Outlier analysis')
    disp('______________________________________________________________')
    
    fprintf('global mean | global mean to var | spatial MAD | ');
else
    verbosestr = 'noverbose';
end

robustz = @(x) (x - mean(x)) ./ mad(x);

gm = mean(dat.dat)';
sm = mad(dat.dat)';

globalmean = abs( robustz( gm )  );
spatialmad = abs( robustz( sm )  );
global_mean_to_var = abs( robustz( gm ./ sm )  );

global_mean_outliers = globalmean > madlim;                       % absolute robust zscore of global mean > limit (default = 3)
global_mean_to_variance_outliers = global_mean_to_var > madlim;   % absolute robust zscore of global mean over Median Abs Deviation > limit
spatialmad_outliers = spatialmad > madlim;                        %spatialmad > mean(spatialmad) + madlim * mad(spatialmad);


% -------------------------------------------------------------------------
% Time series-specific outliers
% -------------------------------------------------------------------------

if dotimeseries
    
    if doverbose
        fprintf('rmssd | ');
    end
    
    % Get RMSSD (DVARS)
    sdiffs = diff(dat.dat')';
    sdiffs = [mean(sdiffs, 2) sdiffs]; % keep in image order
    
    rmssd = ( mean(sdiffs .^ 2) ) .^ .5;    % rmssd - root mean square successive diffs
    rmssd(1) = median(rmssd);               % avoid first time point being very different and influencing distribution and plots.
    rmssd = rmssd';
    
    rmssd_outliers = rmssd > mean(rmssd) + madlim * mad(rmssd);
    
else
    % Omit time series measures -- all false
    
    rmssd = NaN .* zeros(size(global_mean_outliers));
    rmssd_outliers = false(size(global_mean_outliers));
    
end

% -------------------------------------------------------------------------
% Coverage
% -------------------------------------------------------------------------

if doverbose
    fprintf('Missing values | ');
end

desc = descriptives(dat, 'noverbose');
missingvals = desc.images_missing_over_25percent';

if doverbose
    fprintf('%3.0f images \n', sum(missingvals));
    fprintf('\n\n');
end

% -------------------------------------------------------------------------
% Multivariate outliers
% -------------------------------------------------------------------------

% Get mahalanobis outliers
[mahalcov, ~, ~, mahal_cov_outlier_uncorr, mahal_cov_outlier_corr] = mahal(dat, 'noplot', verbosestr);
[mahalcorr, ~, ~, mahal_corr_outlier_uncorr, mahal_corr_outlier_corr] = mahal(dat, 'corr', 'noplot', verbosestr);

if doverbose
    fprintf('Mahalanobis (cov and corr, q<0.05 corrected):\n');
    fprintf('%3.0f images \n', sum(mahal_cov_outlier_corr | mahal_corr_outlier_corr));
end

% -------------------------------------------------------------------------
% Framewise Displacement
% -------------------------------------------------------------------------

if dofd
    % Get FD
    [~, est_outliers, net_mvmt] = framewise_displacement(fd);

    % identify FD outliers based on statistical distribution
    mu = nanmean(net_mvmt);
    sd = nanstd(net_mvmt);
    fd_cdf = cdf('norm',net_mvmt,mu,sd);
    fd_p = 2*min(fd_cdf, 1-fd_cdf);

    fd_corr_out = false(size(fd_p));
    fdr_thr = FDR(fd_p, 0.05);
    if ~isempty(fdr_thr)
        fd_corr_out(fd_p < FDR(fd_p, 0.05) & est_outliers) = true;
    end
    fd_uncorr_out = fd_p < 0.05 | est_outliers;
    
    if doverbose
        fprintf('Framewise Displacement (before and after >.25mm correction):\n');
        fprintf('%3.0f images \n', sum(fd_uncorr_out | fd_corr_out));
    end
end

% -------------------------------------------------------------------------
% Summarize
% -------------------------------------------------------------------------

if dofd
    est_outliers_uncorr = global_mean_outliers | global_mean_to_variance_outliers | rmssd_outliers | spatialmad_outliers | mahal_cov_outlier_uncorr | mahal_corr_outlier_uncorr | fd_corr_out | missingvals;
    est_outliers_corr = global_mean_outliers | global_mean_to_variance_outliers | rmssd_outliers | spatialmad_outliers | mahal_cov_outlier_corr | mahal_corr_outlier_corr | fd_uncorr_out | missingvals;
    
    % Make indicator table
    
    outlier_indicator_table = table(global_mean_outliers, global_mean_to_variance_outliers, missingvals, rmssd_outliers, spatialmad_outliers, mahal_cov_outlier_uncorr, mahal_cov_outlier_corr, mahal_corr_outlier_uncorr, mahal_corr_outlier_corr, fd_uncorr_out, fd_corr_out, est_outliers_uncorr, est_outliers_corr, ...
        'VariableNames', {'global_mean' 'global_mean_to_variance' 'missing_values', 'rmssd_dvars', 'spatial_variability', 'mahal_cov_uncor', 'mahal_cov_corrected', 'mahal_corr_uncor', 'mahal_corr_corrected', 'fd_uncorrected', 'fd_corrected', 'Overall_uncorrected', 'Overall_corrected'});

else
    est_outliers_uncorr = global_mean_outliers | global_mean_to_variance_outliers | rmssd_outliers | spatialmad_outliers | mahal_cov_outlier_uncorr | mahal_corr_outlier_uncorr | missingvals;
    est_outliers_corr = global_mean_outliers | global_mean_to_variance_outliers | rmssd_outliers | spatialmad_outliers | mahal_cov_outlier_corr | mahal_corr_outlier_corr | missingvals;
    
    % Make indicator table
    
    outlier_indicator_table = table(global_mean_outliers, global_mean_to_variance_outliers, missingvals, rmssd_outliers, spatialmad_outliers, mahal_cov_outlier_uncorr, mahal_cov_outlier_corr, mahal_corr_outlier_uncorr, mahal_corr_outlier_corr, est_outliers_uncorr, est_outliers_corr, ...
        'VariableNames', {'global_mean' 'global_mean_to_variance' 'missing_values', 'rmssd_dvars', 'spatial_variability', 'mahal_cov_uncor', 'mahal_cov_corrected', 'mahal_corr_uncor', 'mahal_corr_corrected' 'Overall_uncorrected' 'Overall_corrected'});

end

outliercounts = sum(table2array(outlier_indicator_table))';

summary_table = table(outliercounts, 100 * outliercounts ./ desc.n_images, 'VariableNames', {'Outlier_count' 'Percentage'}, 'RowNames', outlier_indicator_table.Properties.VariableNames);

if doverbose, disp(summary_table); end

% Make score Table

score_table = table(globalmean, global_mean_to_var, spatialmad, rmssd, mahalcov, mahalcorr, ...
    'VariableNames', {'globalmean' 'global_mean_to_var' 'spatialmad' 'rmssd_dvars' 'mahal_cov' 'mahal_corr'});

outlier_tables = struct('outlier_indicator_table', outlier_indicator_table, 'score_table', score_table, 'summary_table', summary_table);

outlier_tables.est_outliers_uncorr = est_outliers_uncorr;
outlier_tables.est_outliers_corr = est_outliers_corr;

% -------------------------------------------------------------------------
% Make indicator matrices
% -------------------------------------------------------------------------

outlier_tables.outlier_regressor_matrix_uncorr = intercept_model(size(dat.dat, 2), find(est_outliers_uncorr));
outlier_tables.outlier_regressor_matrix_uncorr = outlier_tables.outlier_regressor_matrix_uncorr(:, 2:end);  % remove initial intercept column

outlier_tables.outlier_regressor_matrix_corr = intercept_model(size(dat.dat, 2), find(est_outliers_corr));
outlier_tables.outlier_regressor_matrix_corr = outlier_tables.outlier_regressor_matrix_corr(:, 2:end);  % remove initial intercept column


% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
if doplot

    if dobriefplot
        hold on;
        scores = zscore(table2array(score_table));
        maxscore = nanmax(scores, [], 2);
        plot(scores, 'k.-'); %, 'MarkerSize', 4);

        plot(find(est_outliers_uncorr), maxscore(est_outliers_uncorr), '+', 'color', [1 .3 .3], 'MarkerSize', 4, 'LineWidth', 2, 'MarkerFaceColor', [.5 .25 0]);
        plot(find(est_outliers_corr), maxscore(est_outliers_corr), 'ro', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', [1 .5 0]);

        xlabel('Case number');
        ylabel('Scaled outlier criterion scores');

    else % full plot
        create_figure('plot', 3, 2);
        plot(score_table.globalmean); plot(find(global_mean_outliers), score_table.globalmean(find(global_mean_outliers)), 'ro', 'MarkerFaceColor', 'r');
        title('Global mean')

        subplot(3, 2, 2)
        plot(global_mean_to_var); plot(find(global_mean_to_variance_outliers), global_mean_to_var(find(global_mean_to_variance_outliers)), 'ro', 'MarkerFaceColor', 'r');
        title('Global mean to var')

        subplot(3, 2, 3)
        x = spatialmad;
        wh = spatialmad_outliers;
        plot(x); plot(find(wh), x(find(wh)), 'ro', 'MarkerFaceColor', 'r');
        title('Spatialmad')

        subplot(3, 2, 4)
        x = mahalcov;
        wh = mahal_cov_outlier_corr;
        plot(x); plot(find(wh), x(find(wh)), 'ro', 'MarkerFaceColor', 'r');
        title('Mahal cov')

        subplot(3, 2, 5)
        x = rmssd;
        wh = rmssd_outliers;
        plot(x); plot(find(wh), x(find(wh)), 'ro', 'MarkerFaceColor', 'r');
        title('RMSSD')
        xlabel('Case number');


        subplot(3, 2, 6)
        x = mahalcorr;
        wh = mahal_corr_outlier_corr;
        plot(x); plot(find(wh), x(find(wh)), 'ro', 'MarkerFaceColor', 'r');
        title('Mahal corr')
        xlabel('Case number');

    end


end

end % main function


