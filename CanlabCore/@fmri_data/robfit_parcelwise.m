function OUT = robfit_parcelwise(imgs, varargin)
% Robust multiple regression on averages within each parcel in an atlas-class object
%
% - Running parcel-wise saves computation time and reduces multiple comparisons problems
% - Generates t-maps for each predictor in design matrix saved in imgs.X
% - Generates results tables with atlas-labeled regions and montages/surface plots
% - Generates diagnostics on input images, including coverage, inter-subject
% scaling differences, global signal variation across individuals, weights
% across the brain, and more
% - Saves results, t-map object, tables, diagnostic metrics in OUT structure
%
% :Usage:
% ::
%
%      OUT = robfit_parcelwise(imgs, varargin)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2021 Tor Wager
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
%   **imgs:**
%        A dataset of images for 2nd-level analysis
%        Generally 1st-level contrast images, with one image per participant
%
% :Optional Inputs:
%   **names:**
%        followed by cell array of names for each regressor
%
%   **'doplot', [logical flag]:
%        Create plots; default = true. 'noplot' to turn off.
%
%   **'doverbose', [logical flag]:
%        Verbose output; default = true. 'noverbose' to turn off.
%
%   **'csf_wm_covs', [logical flag]:
%        Add global WM and CSF from standard MNI-space masks as covariates
%        default = false 
% 
%   **'remove_outliers', [logical flag]:
%        Remove outliers identified based on mahalanobis distance on either 
%        cov or corr across images at p < 0.05 uncorrected 
%        default = false 
%
% :Outputs:
%
%   **OUT:**
%                  betas: [489×2 double]    Parcels x Predictors regression slopes
%                tscores: [489×2 double]    Parcels x Predictors t-scores
%                pvalues: [489×2 double]    Parcels x Predictors p-values
%              nsubjects: [489×1 double]    Number of subjects analyzed
%                maskvol: [489×1 double]    Vector of parcels run (1) or missing (0)
%                weights: [489×30 double]   Parcels x Subjects robust regression weights
%                    dfe: [489×1 double]    error df for each parcel
%            pthr_FDRq05: [0.0081 0.0128]   P-thresholds for FDR q < .05 for each predictor (map)
%                sig_q05: [489×2 logical]   Parcels x predictors q < 0.05
%         cohens_d_fdr05: [0.4834 0.4458]   Min Cohen's d detectable at FDR q < 0.05
%                  t_obj: [1×1 statistic_image] T-statistic objects (statistic_image), thresholded q < 0.05 FDR initially
%               beta_obj: [1×1 fmri_data]   fmri_data object with regression slopes
%                   mask: [1×1 fmri_data]   Mask of voxels in analyzed parcels
%           resultstable: [2×8 table]       Table object with overall summary of results across the brain
%          group_metrics: [1×1 struct]      QC metrics (see qc_metrics_second_level) 
%     individual_metrics: [1×1 struct]
%          outliers_corr: [30×1 logical]    Mahalanobis outliers p < 0.05 uncorrected
%        outliers_uncorr: [30×1 logical]    Mahalanobis outliers q < 0.05 FDR-corrected
%        ind_quality_dat: [30×6 double]     Summary of some QC metrics
%   region_objects: {[1×X region]  [1×X region]} Region objects containing significant blobs with autolabeled names in .shorttitle
%  contrast_tables_FDR05: {[X×8 table]  [X×8 table]} Table objects with labeled significant regions at q < 0.05 FDR 
%
% :Examples:
% ::
%
% imgs = load_image_set('emotionreg');
%
% Xinfo = importdata(which('Wager_2008_emotionreg_behavioral_data.txt'));
% Xnames = Xinfo.textdata(2);
% imgs.X = zscore(Xinfo.data(:, 2));
%
% % Voxel-wise (non-robust, but robust option is available in regress())
% out = regress(imgs);
% t = threshold(out.t, .05, 'fdr');
%
% % Parcel-wise
% OUT = robfit_parcelwise(imgs);
%
% % Post-run interactive visualization, etc.
% o2 = canlab_results_fmridisplay(get_wh_image(t, 2), 'montagetype', 'full', 'noverbose');
% o2 = removeblobs(o2);
% o2 = addblobs(o2, get_wh_image(t_obj, 1));
%
% Adjust for global intensity and global components and re-run:
% imgs2 = imgs.rescale('l2norm_images');
% OUT = robfit_parcelwise(imgs2);
%
% Run adding global WM and CSF signal covariates:
% OUT = robfit_parcelwise(imgs2, 'csf_wm_covs');
%
% run with global covariates and removing outliers based on Mahalanobis distance:
% OUT = robfit_parcelwise(imgs2, 'csf_wm_covs', 'remove_outliers');
%
% :References:
%   None listed yet.
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    Created by Tor Wager, Feb 2021
% ..

% --------------------------------------------
% Defaults and inputs
% --------------------------------------------
% names = {};
% doverbose = true;
% doplots = true;

ARGS = parse_inputs(varargin{:});

% If you want to distribute arguments back out to variables, use this:

fn = fieldnames(ARGS);

for i = 1:length(fn)
    str = sprintf('%s = ARGS.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% --------------------------------------------
% set up images and parcels
% --------------------------------------------

% Checks for required functions
if isempty(which('mafdr')), error('Sorry, you need the Matlab bioinformatics toolbox on your Matlab path to use the Storey 2002 mafdr function called here.'); end

% Load parcel atlas and initalize object
b = brainpathway('noverbose');

% Make image space match atlas space
imgs = resample_space(imgs, b.region_atlas);

% Trigger extraction and averaging of parcel-wise data
b.voxel_dat = imgs.dat;

datmatrix = double(b.region_dat);
[t, v] = size(datmatrix);


% --------------------------------------------
% Run QC metrics and get image diagnostic info
% --------------------------------------------
if doverbose
         
    dashes = '__________________________________________________________________';
    printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);
    
    printhdr('Input image diagnostic information');
    
end

% mahalanobis distance
[ds, ~, ~, wh_outlier_uncorr, wh_outlier_corr] = mahal(imgs, 'noplot', 'corr');
[dscov, ~, ~, wh_outlier_uncorr_cov, wh_outlier_corr_cov] = mahal(imgs, 'noplot');

[group_metrics, individual_metrics, global_gm_wm_csf_values] = qc_metrics_second_level(imgs); % , values, gwcsf, gwcsfmean, gwcsfl2norm


outliers_uncorr = wh_outlier_uncorr | wh_outlier_uncorr_cov;
outliers_corr = wh_outlier_corr | wh_outlier_corr_cov;


% --------------------------------------------
% set up and check design
% --------------------------------------------
X = imgs.X;

if isempty(X)
    X = ones(t, 1);  % for robust reg
else
    % Make sure intercept is last column
    X = intercept(X, 'end');
end

% Check for non-centered predictors (excluding the intercept)
if any(abs(mean(X)) > (100 * eps) & any(abs(X - mean(X)) > 100 * eps, 1))
    warning('SOME PREDICTORS ARE NOT CENTERED - THE INTERCEPT MAP WILL NOT BE INTERPRETABLE AS THE GROUP MEAN');
    pause(2)
end

k = size(X, 2); % number of maps to estimate - one per regressor

names = cell(1, k);
for i = 1:k-1, names{i} = sprintf('Predictor %d', i); end, names{k} = 'Intercept (Group avg)';

% Replace names if entered
if isfield(ARGS, 'names') && ~isempty(ARGS.names)
    
    names = ARGS.names;
    if length(names) < k
        names{k} = 'Intercept (Group avg)';
    end
    
end


% --------------------------------------------
% Additional covariates and adjustments
% --------------------------------------------

% Add CSF and WM covs if requested

if csf_wm_covs
    
    % Add global WM and CSF to design matrix and names
    X = [X(:, 1:end-1) global_gm_wm_csf_values(:, 2:3) X(:, end)];
    
    names = [names(1:end-1) {'Global_WM' 'Global_CSF'} names(end)];
    
    k = k + 2;
    
    if doverbose
        fprintf('Adjustments: Added Global WM and CSF covariates\n', sum(outliers_uncorr));
    end
    
end


if remove_outliers 
    % Remove outliers (uncorr means uncorrected)
    % We do not pass out imgs, so no need to adjust other fields.
    
    X(outliers_uncorr, :) = [];
    
    datmatrix(outliers_uncorr, :) = [];

    if doverbose
        fprintf('Adjustments: Removed %d outlier images\n', sum(outliers_uncorr));
    end

    t = t - sum(outliers_uncorr);
    
    individual_metrics.csf_to_gm_signal_ratio(outliers_uncorr) = []; 
    individual_metrics.gm_L1norm(outliers_uncorr) = []; 
    individual_metrics.csf_L1norm(outliers_uncorr) = [];
    ds(outliers_uncorr) = []; 
    dscov(outliers_uncorr) = [];
    
    global_gm_wm_csf_values(outliers_uncorr, :) = [];
    
end

% --------------------------------------------
% Plot design matrix
% --------------------------------------------

if doplots && k > 2
    
    create_figure('design matrix', 1, 2);
    
    Xplot = zscore(X);
    Xplot(:, end) = 1;
    
    imagesc(Xplot);
    set(gca, 'XTick', 1:k, 'XTickLabel', format_strings_for_legend(names), 'XTickLabelRotation', 45, 'YDir', 'reverse');
    axis tight
    title('Z-scored Design Matrix');
    
    subplot(1, 2, 2);
    
    plot_correlation_matrix(X(:, 1:end-1), 'names', names(1:end-1), 'nofigure');
    title('Correlations among predictors')
    
    subplot(1, 2, 1);
    colormap(gca, 'summer')
    colormap(gca, 'winter'); colorbar

    drawnow, snapnow
    
end

% --------------------------------------------
% Initialize output arrays
% --------------------------------------------
% Parcels are rows, images are columns

[nsubjects, maskvol, dfe] = deal(zeros(v, 1) .* NaN);     % save numbers of missing values

weights = zeros(v, t) .* NaN;                        % robust reg weights

[betas, tscores, pvalues] = deal(zeros(v, k) .* NaN);

% --------------------------------------------
% Run robust reg for each parcel
% --------------------------------------------

for i = 1:v
    
    % Remove NaNs and run only if we have enough observations
    [wasnan, Xvox, tvox] = nanremove(X, datmatrix(:, i));
    
    isok = sum(~wasnan) > 4;
    
    if isok
        % OK to run
        [bb,stats] = robustfit(Xvox, tvox, 'bisquare', [], 'off');
        
        w = naninsert(wasnan, stats.w);
        
    else
        % empty voxel: not enough data
        disp('THIS CODE SHOULD NEVER RUN B/C EXCLUDING BAD VOXELS ABOVE');
        bb = NaN .* ones(k, 1);
        stats = struct('t', bb, 'p', ones(k, 1));
        w = NaN .* ones(N, 1);
        
    end
    
    % Redistribute to maps
    betas(i, :) = bb';
    tscores(i, :) = stats.t';
    pvalues(i, :) = stats.p';
    
    dfe(i, 1) = stats.dfe;
    nsubjects(i, 1) = sum(~wasnan);
    maskvol(i, 1) = isok;
    
    weights(i, :) = w;
    
end % loop through nodes

OUT = struct('betas', betas, 'tscores', tscores, 'pvalues', pvalues, 'nsubjects', nsubjects, 'maskvol', maskvol, 'weights', weights, 'dfe', dfe, 'datmatrix', datmatrix); % @lukasvo76 added datmatrix to OUT struct for flexible plotting

%% --------------------------------------------
% FDR correction
% --------------------------------------------

FDRq = zeros(v, k);
pthr = zeros(1, k);

for i = 1:k
    % for each map
    [FDRq(:, i), ~, pIO] = mafdr(pvalues(:, i));
    
    if pIO > .99
        % all P-values 1, mafdr may not be suitable
        disp('Warning: Prior prob of sig P-values is near 1. Using B-H FDR');
        pthr_i = FDR(pvalues(:, i), .05);
        
    else
        % use mafdr
        % P-threshold for FDR q < 0.05 for each map
        pthr_i = max(pvalues(FDRq(:, i) < 0.05, i));
    end
    
    if isempty(pthr_i)
        pthr(i) = Inf;
    else
        pthr(i) = pthr_i;
    end
end

sig_q05 = FDRq < 0.05;

% Other interesting metrics to save

OUT.pthr_FDRq05 = pthr;
OUT.pthr_FDRq05_descrip = 'Highest P-value significant at FDR q < 0.05; threshold at p <= pthr';
OUT.sig_q05 = sig_q05;

maxn = max(OUT.nsubjects);
OUT.cohens_d_fdr05 = tinv(1 - OUT.pthr_FDRq05, maxn - k) ./ sqrt(maxn - k);
OUT.cohens_d_fdr05_descrip = 'Min Cohen''s d detectable at FDR q < 0.05';

% QC metrics
OUT.group_metrics = group_metrics;
OUT.individual_metrics = individual_metrics;
OUT.global_gm_wm_csf_values = global_gm_wm_csf_values;
OUT.outliers_corr = outliers_corr;
OUT.outliers_uncorr = outliers_uncorr;

Global_Weight_Mean = mean(OUT.weights)';
global_GM = global_gm_wm_csf_values(:, 1);
global_WM = global_gm_wm_csf_values(:, 2);
global_CSF = global_gm_wm_csf_values(:, 3);
CSF_to_GM_ratio = OUT.individual_metrics.csf_to_gm_signal_ratio';
GM_L1_norm = OUT.individual_metrics.gm_L1norm;
CSF_L1_norm = OUT.individual_metrics.csf_L1norm;
Mahal_corr = ds;
Mahal_cov = dscov;

OUT.ind_quality_dat = table(Global_Weight_Mean, global_GM, global_WM, global_CSF, GM_L1_norm, CSF_L1_norm, CSF_to_GM_ratio, Mahal_corr, Mahal_cov);

% --------------------------------------------
% Transform back to voxelwise output maps
% --------------------------------------------

OUT.t_obj = parcel_stats2statistic_image(b.region_atlas, tscores, pvalues, dfe, sig_q05);

OUT.beta_obj = parcel_data2fmri_data(b.region_atlas, betas);

% --------------------------------------------
% Create mask and nsubjects
% --------------------------------------------

mask = parcel_data2fmri_data(b.region_atlas, ones(v, 1));
OUT.mask = mask;

OUT.nsubjects_obj = parcel_data2fmri_data(b.region_atlas, nsubjects);

% mask = get_wh_image(OUT.t_obj, 1);
% mask = threshold(mask, 1 - eps, 'unc', 'noverbose');
% mask = remove_empty(mask);
% mask.dat = ones(size(mask.dat));
% mask.p = .001 * ones(size(mask.dat)); % kludge to avoid region() using p-values

% --------------------------------------------
% Print report & summary results table
% --------------------------------------------

resultstable = table;
resultstable.Properties.Description = 'Parcel-wise robust regression';
resultstable.maxT = max(OUT.tscores)';
resultstable.minP = min(OUT.pvalues)';

resultstable.sig05 = sum(OUT.pvalues < 0.05)';
resultstable.sig005 = sum(OUT.pvalues < 0.005)';
resultstable.sig001 = sum(OUT.pvalues < 0.001)';
resultstable.sigFDR05 = sum(OUT.sig_q05)';
resultstable.p_thr_FDR05 = OUT.pthr_FDRq05';
resultstable.min_d_FDR05 = OUT.cohens_d_fdr05';

resultstable.Properties.RowNames = names;
OUT.resultstable = resultstable;

if doverbose
    
    % Print to screen

    printhdr(resultstable.Properties.Description);
    disp(resultstable);
    disp(' ')
    disp('sig*: Significant parcels at given threshold (p < 0.05 two-tailed, q < 0.05 FDR, etc.)');
    disp('p_thr_FDR05: P-value threshold to achieve q < 0.05 FDR-corrected for each predictor');
    fprintf('min_d_FDR05: %s', OUT.cohens_d_fdr05_descrip);
    disp('dashes')
    disp(' ')
    
end


%%
% --------------------------------------------
% Print montages and tables of regions at q < 0.05 FDR
% --------------------------------------------
if doverbose
    
    printhdr('Tables of regions at q < 0.05 FDR');
    
    OUT.region_objects = cell(1, k);
    
end

for i = 1:k
    
    if doverbose
        printhdr(sprintf('Predictor %d: %s', i, names{i}));
    end
    
    if doplots && ~isinf(OUT.pthr_FDRq05(i)) % if we have some results to show
        
        canlab_results_fmridisplay(get_wh_image(OUT.t_obj, i), 'montagetype', 'full', 'noverbose');
        
        set(gcf, 'Name', names{i}, 'NumberTitle', 'off');
        
    end
    
    if doverbose
        
        r = region(get_wh_image(OUT.t_obj, i), 'noverbose');
        
        if isempty(r)
            OUT.contrast_tables_FDR05{i} = 'No significant results at FDR q < 0.05';
            disp(OUT.contrast_tables_FDR05{i})
            OUT.region_objects{i} = region();
            
        else
            [posr, negr, OUT.contrast_tables_FDR05{i}] = table(r, 'nolegend');
            OUT.region_objects{i} = [posr negr];
        end
        
        disp(' ')
        
        
        
    end

end


% --------------------------------------------
% Additional diagnostic info
% --------------------------------------------
    
if doplots
    
    % --------------------------------------------
    % Mask Figure
    % --------------------------------------------
    
    create_figure('mask'); axis off; montage(OUT.mask, 'color', [0 .5 0], 'trans', 'noverbose');
    
    drawnow, snapnow;
    
    % --------------------------------------------
    % Data Figure
    % --------------------------------------------
    plot(imgs);
    
    drawnow, snapnow;
    
    % --------------------------------------------
    % Weights and diagnostics figure
    % --------------------------------------------
    
    create_figure('weights and metrics', 2, 2);
    xlabel('Image'); ylabel('Weights');
    errorbar(mean(OUT.weights), std(OUT.weights), 'bo', 'MarkerFaceColor', [0 0 .5])
    title('Mean weights across parcels (s.d. error bars) per image');
    axis tight; 
    
    subplot(2, 2, 2);
    imagesc(OUT.weights);
    xlabel('Image'); ylabel('Parcel');
    title('Weights by parcel');
    colorbar;
    axis tight; set(gca, 'YDir', 'Reverse');
    
    
    subplot(2, 2, 3);
    xlabel('Image'); ylabel('Z(Weights)');
    errorbar(zscore(mean(OUT.weights)), ste(OUT.weights), 'bo-', 'MarkerFaceColor', [0 0 .5], 'LineWidth', 2)
    title('Mean weights (s.e. error bars) and quality metrics');
    plot(zscore(OUT.individual_metrics.gm_L1norm), 'LineWidth', 2);
    plot(zscore(OUT.individual_metrics.csf_L1norm), 'LineWidth', 2);
    plot(zscore(ds), 'LineWidth', 2);
    plot(zscore(dscov), 'LineWidth', 2);
    legend({'Z(Weights)' 'Z(GM L1 norm)' 'Z(CSF L1 norm)' 'Mahal corr dist' 'Mahal cov dist'});
    axis tight; 
    
    % mark off who are outliers
    wh_out = find(OUT.outliers_uncorr);
    for i = 1:length(wh_out)
        
        hh = plot_vertical_line(wh_out(i));
        set(hh, 'Color', 'r', 'LineStyle', '--');
        
        if i == 1
                legend({'Z(Weights)' 'Z(GM L1 norm)' 'Z(CSF L1 norm)' 'Mahal corr dist' 'Mahal cov dist' 'Mah. outliers p<.05 uncor'});
        end
    end
    
    subplot(2, 2, 4)
    plot_correlation_matrix(datmatrix, 'dofigure', false);
    title('inter-parcel correlations across images');
    drawnow, snapnow;
    
    
end % doplots

end % function




function ARGS = parse_inputs(varargin)

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_scalar = @(x) validateattributes(x, {'numeric' 'logical'}, {'nonempty', 'scalar'});

% valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector
%
% % Validation: Region object, structure, or [x1 x2 x3] triplet
% valfcn_custom = @(x) isstruct(x) || isa(x, 'region') || (~isempty(x) && all(size(x) - [1 3] == 0) && all(isnumeric(x)));
%
% % Validation: [x1 x2 x3] triplet
% valfcn_xyz = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'size', [1 3]});
%
% valfcn_logical = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', '>=', 0, '<=', 1}); % could enter numeric 0,1 or logical


% Required inputs
% ----------------------------------------------------------------------
% p.addRequired('x', valfcn_custom);
% p.addRequired('y', valfcn_custom);

% Optional inputs
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

% p.addParameter('color', [.9 .2 0], valfcn_xyz);
p.addParameter('names', {}, @iscell); % can be scalar or vector
p.addParameter('doverbose', true, valfcn_scalar);
p.addParameter('doplots', true, valfcn_scalar);
p.addParameter('csf_wm_covs', false, valfcn_scalar);
p.addParameter('remove_outliers', false, valfcn_scalar);



% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);

p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs
