function obj_denoised = denoise_timeseries_pipeline(obj, TR, HP_cutoff_sec, mvmtfname, varargin)
% denoise_timeseries_pipeline - Prepare 4-D time series images for connectivity or model-based analysis.
%
% :Usage:
% ::
%     obj_denoised = denoise_timeseries_pipeline(obj, TR, HP_cutoff_sec, mvmtfname, [optional inputs])
%
% :Inputs:
%
%   **obj:**
%        An fmri_data object containing a 4-D time series image.
%
%   **TR:**
%        Scalar. The repetition time (in seconds) of the acquisition.
%
%   **HP_cutoff_sec:**
%        Scalar. High-pass filter cutoff in seconds.
%
%   **mvmtfname:**
%        String. Full path to the movement parameter file (rp_.txt) produced after realignment.
%
% :Optional Inputs:
%
%   **'plot':** [logical]
%        If true (default = true), display intermediate figures.
%
%   **'verbose':** [logical]
%        If true (default = true), print detailed processing information.
%
%   **'save':** [logical]
%        If true (default = false), save denoised object in the original file location in obj.fullpath, with '_denoised.mat' extension.
%
%   **'pca_denoise':** [logical]
%        If true (default = false), perform additional PCA-based denoising.
%
%   **'images_per_run':** [numeric vector]
%        A vector indicating the number of image volumes per run. If provided,
%        indicator regressors for each run will be generated and added to the covariates.
%       Note: If your object already contains this info in
%       .images_per_session, this info will be used.
%
% :Outputs:
%
%   **obj_denoised:**
%        An fmri_data object containing the denoised 4-D time series image.
%        The denoised object includes:
%           - .dat: The residual data after regressing out nuisance covariates.
%           - .covariates: The matrix of nuisance regressors used in the denoising.
%           - .metadata_table: A table summarizing the nuisance regressors, including movement regressors,
%                              run intercepts (if applicable), spike/outlier regressors, CSF signal,
%                              and high-pass filter components.
%
% :Steps:
%
%   1. Load the movement parameters file (mvmtfname) and compute 24 movement-related regressors.
%   2. If 'images_per_run' is provided, create run-specific intercept regressors and add these to the covariate matrix.
%   3. Identify spike and high framewise-displacement outlier images using outliers(), and append the resulting regressors.
%   4. Extract the mean CSF signal using extract_gray_white_csf() and add it as a regressor.
%   5. (Optional) Perform PCA-based denoising if 'pca_denoise' is true.
%   6. Apply a high-pass filter to the time series using use_spm_filter with TR and HP_cutoff_sec.
%   7. Rescale the combined nuisance covariates for visualization.
%   8. Visualize the final covariate matrix if requested
%   9. Regress out the nuisance covariates from the 4-D time series (using regress with residual and grand mean scaling).
%  10. Attach the covariates and metadata table to the denoised object.
%  11. Plot the denoised object data (if requested)
%  12. Save the denoised object to disk if requested.
%
% :Notes on processes:
% - Uses the 'residual' option in fmri_data.regress() to return residuals after removing covariates
% - Use the 'add_voxelwise_intercept' option to preserve the mean signal in
% each voxel. This doesn't matter if we are only calculating correlations,
% but can matter if we are doing anything that requires or interacts with
% maps of the overall level of signal in each voxel, like spatial
% correlations or SPM analyses with implicit masking.
% - use the 'grandmeanscale' option to rescale the entire 4-D dataset
% (usually run) to a grand mean of 100.  This doesn't matter if we're only
% calculating within-run correlations, but matters if we are going to stack
% data across runs that may be on a different scale.
% Scale run-wise, using obj.images_per_session to rescale each run to a grand mean of 100 
% regress() will use obj.rescale('session_grand_mean_scaling_spm_style') to scale each run
%
% :Examples:
% ::
%     % Load the key 4-D image file into an object:
%     fname = which('swrsub-sid001567_task-pinel_acq-s1p2_run-03_bold.nii.gz');
%     obj = fmri_data(fname);
%
%     % Get the TR from the JSON file:
%     json_struct = jsondecode(fileread(which('sub-sid001567_task-pinel_acq-s1p2_run-03_bold.json')));
%     TR = json_struct.RepetitionTime;
%
%     % Specify the movement parameter file:
%     mvmtfname = which('rp_sub-sid001567_task-pinel_acq-s1p2_run-03_bold.txt');
%
%     % Run the denoising pipeline:
%     obj_denoised = obj.denoise_timeseries_pipeline(TR, 128, mvmtfname, 'plot', true, 'verbose', true);
%
% :See also:
%   movement_regressors, intercept_model, outliers, extract_gray_white_csf, pca, use_spm_filter, regress
%

% Author: Tor Wager
% Date: 2025-Mar-27
% License: GNU General Public License v3 or later

%% Parse Optional Inputs
p = inputParser;
addParameter(p, 'plot', true, @(x) islogical(x));
addParameter(p, 'verbose', true, @(x) islogical(x));
addParameter(p, 'save', false, @(x) islogical(x));
addParameter(p, 'pca_denoise', false, @(x) islogical(x));
addParameter(p, 'images_per_run', [], @(x) isnumeric(x));
parse(p, varargin{:});
doPlot = p.Results.plot;
verbose = p.Results.verbose;
doSave = p.Results.save;
pca_denoise = p.Results.pca_denoise;
images_per_run = p.Results.images_per_run;

% Scale run-wise, using obj.images_per_session to rescale each run to a grand mean of 100 
% regress() will use this to scale each run  
if ~isempty(images_per_run)
    obj.images_per_session = images_per_run;
end

    

if verbose, disp('Starting denoise_timeseries_pipeline'); end

%% Step 1: Load Movement Regressors
if verbose, disp('Loading movement regressors...'); end
[mvmt_matrix, mvmt_regs_24, mvmt_table] = movement_regressors(mvmtfname);
covs = mvmt_regs_24;

%% Step 2: Add Run-Specific Regressors (if applicable)

if ~isempty(obj.images_per_session)

    total_images = size(obj.dat, 2);
    if sum(obj.images_per_session) ~= total_images
        error('The sum of images_per_run (%d) does not equal the total number of images (%d).', sum(obj.images_per_session), total_images);
    end

    % Create run intercept regressors using intercept_model
    run_indicators = intercept_model(obj.images_per_session);
    covs = [covs run_indicators];

    % Create a table with variable names 'Run1', 'Run2', ..., and add to cov_table.
    runNames = arrayfun(@(i) sprintf('Run%d', i), 1:length(obj.images_per_session), 'UniformOutput', false);
    run_table = array2table(run_indicators, 'VariableNames', runNames);
    cov_table = [mvmt_table, run_table];

else
    cov_table = mvmt_table;
end

%% Step 3: Add Outlier Regressors

if verbose, disp('Computing outlier regressors...'); end

[est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(obj, 'fd', mvmt_matrix);

if verbose, fprintf('Number of outlier volumes (uncorrected): %d\n', sum(est_outliers_uncorr)); end

covs = [covs outlier_tables.outlier_regressor_matrix_uncorr];

out_table = array2table(outlier_tables.outlier_regressor_matrix_uncorr);
out_table.Properties.VariableNames = strrep(out_table.Properties.VariableNames, 'Var', 'Out');

cov_table = [cov_table out_table];

%% Step 4: Add Mean CSF Regressor

if verbose, disp('Extracting mean CSF signal...'); end

[gwcsfvalues, gwcsfcomponents] = extract_gray_white_csf(obj);
meancsf = gwcsfvalues(:, 3);
covs = [covs meancsf];
cov_table = addvars(cov_table, meancsf, 'NewVariableNames', 'MeanCSF');

%% Step 5: Optional PCA Denoising

if pca_denoise
    if verbose, disp('Performing PCA denoising...'); end
    [component_scores, eigenmap_obj, explained] = pca(obj, 'k', 20, 'noplot');

    r2 = regress_component_scores(component_scores, mvmt_regs_24);

    wh = r2 > .7;

    if verbose, fprintf('Number of components with r^2 > 0.7 removed: %d\n', sum(wh)); end

    keep = ~wh;  % logical index of components to keep

    X_recon = (obj.dat' * eigenmap_obj.dat(:, keep)) * eigenmap_obj.dat(:, keep)';  % Reconstruct X without the excluded components
    obj.dat = X_recon';

end

%% Step 6: High-Pass Filter

if verbose, disp('Applying high-pass filter...'); end

[~, ~, KH] = use_spm_filter(TR, size(obj.dat, 2), 'none', 'specify', HP_cutoff_sec);
covs = [covs KH];

KH_table = array2table(KH);
KH_table.Properties.VariableNames = strrep(KH_table.Properties.VariableNames, 'Var', 'HPfilt');
cov_table = [cov_table KH_table];

%% Step 7: Rescale Covariates

covs = covs ./ max(abs(covs), [], 1);

%% Step 8: Visualize Covariates (if requested)

if doPlot

    create_figure('all covs');
    imagesc(covs);
    set(gca, 'YDir', 'Reverse'); 
    axis tight;
    colormap(colormap_tor([0 0 1], [1 1 0], [.2 .5 1], [.5 .5 .5], [1 .5 .2]));
    title('Nuisance covariates removed')

end

%% Step 9: Regress Out Nuisance Covariates

if verbose, disp('Regressing out nuisance covariates...'); end

obj.X = covs;

% Scale run-wise, using obj.images_per_session to rescale each run to a grand mean of 100 
% regress() will use obj.rescale('session_grand_mean_scaling_spm_style') to scale each run 

regress_output = regress(obj, 'residual', 'grandmeanscale', 'add_voxelwise_intercept', 'noverbose');
obj_denoised = regress_output.resid;

%% Step 10: Attach Covariates and Metadata

obj_denoised.covariates = covs;
obj_denoised.metadata_table = cov_table;

if verbose
    disp('Denoised metadata table:');
    disp(obj_denoised.metadata_table);
end

%% Step 11: Plot the denoised object data (if requested)

if doPlot

    plot(obj_denoised)

end

%% Step 12: Save the Denoised Object (if requested)

obj_denoised = enforce_variable_types(obj_denoised);

if doSave
    [dd, ff, ~] = fileparts(obj_denoised.fullpath(1, :));
    [~, fn] = fileparts(ff);
    fn = [fn '_denoised.mat'];
    save(fullfile(dd, fn), 'obj_denoised');

    if verbose
        fprintf('Denoised object saved as: %s\n', fullfile(dd, fn));
    end
end


if verbose, disp('denoise_timeseries_pipeline finished.'); end

end % Main function





function r2 = regress_component_scores(component_scores, mvmt_regs_24)
% regress_component_scores Regress each column of component_scores on mvmt_regs_24.
%
% :Usage:
% ::
%     r2 = regress_component_scores(component_scores, mvmt_regs_24)
%
% :Inputs:
%
%   **component_scores:** [n x k] numeric matrix.
%        Each column contains a set of component scores (e.g., from PCA).
%
%   **mvmt_regs_24:** [n x v] numeric matrix.
%        Matrix of movement regressors (e.g., 24 motion parameters).
%
% :Outputs:
%
%   **r2:** [1 x k] numeric vector.
%        The coefficient of determination (R^2), representing the variance 
%        explained by mvmt_regs_24 (with intercept) in each column of component_scores.
%
% :Examples:
% ::
%     % Example: Compute R^2 for each PCA component regressed on movement regressors.
%     n = 100;
%     k = 5;
%     v = 6;
%     component_scores = randn(n, k);
%     mvmt_regs_24 = randn(n, v);
%     r2 = regress_component_scores(component_scores, mvmt_regs_24);


    % Check that the number of rows in both matrices is equal.
    if size(component_scores, 1) ~= size(mvmt_regs_24, 1)
        error('The number of rows in component_scores must equal the number of rows in mvmt_regs_24.');
    end
    
    n = size(component_scores, 1);
    k = size(component_scores, 2);
    
    % Create design matrix with intercept.
    X = [ones(n, 1) mvmt_regs_24];
    
    r2 = zeros(1, k);
    for i = 1:k
        y = component_scores(:, i);
        % Compute regression coefficients.
        b = regress(y, X);
        yhat = X * b;
        % Compute sum of squared errors and total sum of squares.
        SSE = sum((y - yhat).^2);
        SST = sum((y - mean(y)).^2);
        if SST == 0
            r2(i) = NaN;
        else
            r2(i) = 1 - SSE / SST;
        end
    end
end