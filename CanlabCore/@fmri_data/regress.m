function out = regress(obj, varargin)
% Multiple regression with an fmri_data object (obj), probj.objedicting brain data with a design matrix stored in obj.X (or vice versa)
%
% Regress obj.X on obj.dat at each voxel, and return voxel-wise statistic
% images. Each column of obj.X is a predictor in a multiple regression,
% and the intercept is the last column. Intercept will automatically be
% added if not detected unless 'nointercept' is specified.
%
% Key pointers - Using output:
% - Output structure (regression_results_ols) contains beta and t images in statistic_image objects
% - If contrasts are entered, output structure contains contrast images and contrast t images in statistic_image objects
% - T images are thresholded, but beta/contrast images are not. These can be re-thresholded with the threshold( ) method
% - statistic_image objects can be visualized with any CANlab display methods, e.g., surface( ), montage( ), and orthviews( )
% - Use the get_wh_image( ) method to select specific beta, t, contrast, contrast t images from statistic_image objects
% - Enter condition names and contrast names as cell vectors of strings Labels are saved in the .image_labels property in output objects
% - Some diagnostics, including VIFs and leverages, are automatically checked and returned
%
% Key pointers - Input options:
% - There is an option for robust regression, specified using the 'robust' flag
% - There is an option for autoregressive (AR) time series modeling, specified using the 'AR' keyword
% - There is an option for voxel-varying covariates, e.g., post-treatment controlling for pre-treatment in same voxel. 'robust' mode only - not in OLS
% - For first-level time series models (only), the grandmeanscale option is recommended to increase homogeneity in scale across participants
% - The regress( ) method does not use covariates field of fmri_data(). You must include covariates manually in obj.X.
%
% - The regress( ) method can also create a map of brain regions that predict the obj.Y vector using the 'brainony' option.
%   This is essentially a univariate version of the 'predict' command.  Warning: this is very slow as it loops
%   through all voxels.
%
% :Usage:
% ::
%
%    regression_results = regress(obj, varargin)
%
% :Inputs:
%  **obj:**
%        should be an fmri_data object with X field defined.
%        obj.X can be a design_matrix() object.
%
% :Optional Inputs:
%  **[threshold, 'unc']:**
%        p-value threshold string indicating threshold type
%        (see help statistic_image.threshold for options)
%
%  **robust:**
%        Run a robust regression (default is OLS).  Robust is considerably
%        slower than OLS
%
%  **grandmeanscale:**
%        Scale overall grand mean to a value of 100. Intended to reduce inter-subject
%        variability in 1st-level analysis (single-subject) when doing
%        multi-subject group analyses on resulting contrasts/beta images.
%        Assumes mask and overall brain size are consistent across replicates (e.g., participants, in a first-level analysis)
%        Do not use for 2nd-level or standard single-level analyses.
%
%  **C**
%       Followed by a contrast matrix, each column is a contrast across conditions/events
%       [k x c] matrix, where k = size(X, 2) and c is number of contrasts
%
%       Note: k must be the number including the intercept. If your input
%       model X does not include an intercept, contrast matrix must have
%       one more element than input X has rows, to account for the added
%       intercept.
%
%  **nointercept:**
%        Do not add intercept to model
%
%  **display, display_results:**
%        Show thresholded results usin orthviews (if < 10 regressors)
%
%  **nodisplay:**
%        Do not plot thresholded results using orthviews
%
%  **brainony:**
%        univariate approach to predict obj.Y from brain data
%
%  **residual:**
%        Output residual as fmri_data() object
%
%  **noverbose:**
%        Suppress verbose outputs
%
%  **variable_names:** (or 'names')
%       Followed by a cell array of variable/regressor names, for
%       non-intercept regressors
%
%  **contrast_names:**
%       Followed by a cell array of contrast names
%
%  **analysis_name:**
%       Followed by a string with a name/description for this analysis.
%
%  **covdat:**
%       Followed by an fmri_data object, or a cell array of fmri_data
%       objects. These serve as voxel-varying covariates -- each voxel can have
%       a unique set of covariates. One common use case might be to regress
%       post-treatment data on group + pre-treatment data. NOTE: this field
%       only works in 'robust' mode, and is not implemented for OLS regression.
%
% :Outputs:
%
%  **regression_results:**
%        A structure containing stats_img and fmri_data objects.
%        In addition to the main outputs below, the regression_results structure also has
%        fields for input_parameters, the design matrix (X), variable
%        names, and warnings.
%
%  **regression_results.b:**
%        stats_img object of beta values estimated from regression
%
%  **regression_results.t:**
%        stats_img object of t-values with input threshold
%
%  **regression_results.df:**
%        fmri_data object of degrees of freedom
%
%  **regression_results.sigma:**
%        fmri_data object of variance of residual
%
%  **regression_results.residual:**
%        fmri_data object of residual data after model has been regressed out (optional).
%
%  **regression_results.diagnostics:***
%        A structure containing VIFs and leverage values for the design matrix
%        regression_results.diagnostics.Variance_inflation_factors = VIFs
%        regression_results.diagnostics.Leverages = leverage values
%
% :Examples:
% ::
%
% % ---------------------------------------------------------------------
% % A simple 2nd-level group analysis
% % ---------------------------------------------------------------------
% % 1. Load a sample dataset with 30 participants from an emotion regulation
% %    task (Wager et al. 2008, Neuron)
% obj = load_image_set('emotionreg');
%
% % 2. Define predictor matrix from metadata_table
% %    Here this is a single column.
% %    You can add an intercept, but if you do not, it will be added automatically
% %    If you mean-center the regressor(s) in X, you will be able to interpret
% %    the intercept map as the activation for the average image (e.g.,
% %    participant in a 2nd-level analysis).
% obj.X = obj.metadata_table.Reappraisal_Success;
% obj.X = obj.X - mean(obj.X);                      % mean-center
% obj.X(:, end + 1) = 1;                            % intercept
%
% % 3. Run the regression and create a plot of the results
% regression_results = regress(obj, .05, 'unc');
%
% %    Output includes statistic_image objects with t maps. These contain
% %    t-values and p-values for each voxel, and a list of which voxels are
% %    significant.
% regression_results.t
%
% % 4. Naming: It is helpful to name the regressors
% % If you do, the output statistic_object maps will be labeled
% regression_results = regress(obj, .05, 'unc', 'names', {'Reapp_Success' 'Intercept'}, 'analysis_name', 'Emotion Regulation');
% regression_results.t.image_labels
%
% % 5. Visualizing thresholded t-maps
% %    There is one image per regressor in X, including the intercept. You can
% %    select, visualize, and re-threshold these:
% montage(regression_results.t)
% t_obj_reapp_success = get_wh_image(regression_results.t, 1);
% create_figure('surface'); surface(t_obj_reapp_success);
%
% t_obj_reapp_success = threshold(t_obj_reapp_success, .005, 'unc');
% create_figure('surface at P < 0.005'); surface(t_obj_reapp_success);
% create_figure('montage at P < 0.005'); axis off; montage(t_obj_reapp_success)
%
% % 6. Add a 2nd-level nuisance covariate
% % It is often helpful to add covariates related to artifacts or effects
% of no interest. In this case we will extract and add mean CSF as a
% covariate. We find that mean gray matter, white matter, and CSF are all
% highly correlated. But there is likely no real signal in CSF, so much of
% what it captures likely reflects image-wide confounds and noise.
%
% gwcsf = obj.extract_gray_white_csf;
% obj.X(:, 2) = gwcsf(:, 3);  % add mean CSF as a covariate
% obj.X(:, 3) = 1;
% regression_results = regress(obj, 'names', {'Reapp_Success' 'CSF_mean' 'Intercept'}, 'analysis_name', 'Emotion Regulation');
% montage(regression_results.t)
%
% % 7. Apply a gray-matter mask before thresholding
% % When correcting for multiple comparisons, it is often advantageous to
% % correct within gray matter voxels, as these are ones we hypothesize. This
% % can make results more significant.
% obj = apply_mask(obj, fmri_data(which('gray_matter_mask.nii')));
% regression_results.t = threshold(regression_results.t, 0.05, 'fdr');
% montage(regression_results.t)
%
% % 8. Look for outliers in several ways
% We are interested in this example in mahal_corr, outlier images that do not
% correlate with the other images.
% We run the outliers in 'notimeseries' mode, which omits some time
% series-specific calculations as this is a 2nd-level analysis
% [~, wh_outliers, outlier_table] = obj.outliers('notimeseries');
% find(wh_outliers)  % these cases are outliers
% %
% % Now we can exclude these and re-run
% obj_cleaned = get_wh_image(obj, ~wh_outliers);
% regression_results = regress(obj_cleaned, 'names', {'Reapp_Success' 'CSF_mean' 'Intercept'}, 'analysis_name', 'Emotion Regulation');
% regression_results.t = apply_mask(regression_results.t, fmri_data(which('gray_matter_mask.nii')));
% regression_results.t = threshold(regression_results.t, 0.01, 'unc');
% montage(regression_results.t)
%
% ---------------------------------------------------------------------
% Other examples:
%
%    % Run regression with liberal threshold
%    regression_results = regress(obj, .05, 'unc');
%
%    % Run regression with conservative threshold and save residual
%    regression_results = regress(obj, .001, 'unc', 'residual);
%
%    % Run robust regression with fdr threshold
%    regression_results = regress(obj, .05, 'fdr','robust');
%
%    % Run a regression predicting behavior from brain at liberal threshold
%    regression_results  = regress(data_comb, .05, 'unc', 'brainony')
%
%    % Re-threshold at different values
%    regression_results.t = threshold(regression_results.t, .05, 'fdr');
%    regression_results.t = threshold(regression_results.t, .001, 'unc');
%
%    % Re-display results of thresholding
%    orthviews(regression_results.t);
%
%    % Write out beta image(s) to current directory
%    regression_results.b.fullpath = fullfile(pwd,'beta.nii');
%    write(regression_results.b)
%
%    % Plot diagnostics
%    figure; subplot(1,2,1); title('VIFs')
%    plot(regression_results.diagnostics.Variance_inflation_factors);
%
%    subplot(1,2,2); title('Leverage of each observation')
%    plot(regression_results.diagnostics.Leverages);
%
%   % Run with options:
%   regression_results = regress(obj, 'variable_names', names, 'analysis_name', 'Pinel localizer 1st-level GLM', 'noverbose');
%
%   regression_results_ols = regress(obj, 'variable_names', names, ...
%     'analysis_name', 'Pinel localizer 1st-level GLM', ...
%     'C', C.weights, 'contrast_names', C.names);
%
% Show some results for selected individual conditions
% % Thresholded t map for the first condition
% t = get_wh_image(regression_results_ols.t, 1);
%
% create_figure('montage'); axis off
% display_obj = montage(t);
%
% figure; surface(t);
%
% % Montages for thresholded t images, condition [2 4 6 8]
%
% t = get_wh_image(regression_results_ols.t, [2 4 6 8]);
%
% create_figure('montage2'); axis off
% display_obj = montage(t);
% Show  results for contrasts
% t = get_wh_image(regression_results_ols.con_t, 1:3);
%
% create_figure('montage3'); axis off
% display_obj = montage(t);
%
% t = get_wh_image(regression_results_ols.con_t, 4:6);
%
% create_figure('montage4'); axis off
% display_obj = montage(t);

% ..
%    Copyright (c) 2015 Tor Wager & Luke Chang
%
%    Permission is hereby granted, free of charge, to any person obtaining a
%    copy of this software and associated documentation files (the "Software"),
%    to deal in the Software without restriction, including without limitation
%    the rights to use, copy, modify, merge, publish, distribute, sublicense,
%    and/or sell copies of the Software, and to permit persons to whom the
%    Software is furnished to do so, subject to the following conditions:
%
%    The above copyright notice and this permission notice shall be included
%    in all copies or substantial portions of the Software.
%
%    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
%    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
%    DEALINGS IN THE SOFTWARE.
% ..

% ..
%    Programmers' Notes:
%    c Tor Wager, Dec 2010
%    Edited by Luke Chang, 9/27/2012 to add optional input to reverse X & Y (i.e., create a map of voxels that predict the behavioral variable)
%    Edited by Luke Chang, 9/28/2012 to add optional input to run robust regression for brainony
%    Edited by Luke Chang, 10/24/2012 to save residuals (i.e., regression_results.r), which is helpful for denoising an image
%    Edited by Luke Chang, 3/26/2013 to add optional input to not add an intercept - allows for more flexible modeling options
%    Code completely refactored by Luke Chang 2/24/25
%    Verbose option updated by Tor, 7/2015
%    Help, some outputs, robust % update: added by Tor, 5/2021
%   Edited by Yoni Ashar - Voxel-varying covariates, @2023-2024
%   Edited by Ke Bo 2/25/2025 to add an option input to include fit general least square regression with AR model. AR order can be defined by users.
% . Example:  Out = regress(fmri_data,'AR',3);
%   Edited by Tor Wager, Mar 2025, to add 'add_voxelwise_intercept' option for residuals
%   and update help with more examples

% ..
%    ---------------------------------------------------------------------
%    Defaults
%    ---------------------------------------------------------------------
% ..
inputargs = {.001, 'uncorrected'}; % default options for thresholding
do_display = false;
brain_is_outcome = 1; %else do brain on Y
do_robust = 0;
do_intercept = 1;
do_resid = 0;
doverbose = true;
variable_names = {};
contrast_names = {};
analysis_name = '';
C = [];
grandmeanscale = false;
covdat = [];   % covariates for each vox. initialize to empty (no voxel-varying covariates)

do_AR = false;
ar_order = 0;  % default AR order (0 means not used)
add_voxelwise_intercept = false;

% ---------------------------------------------------------------------
% Parse Inputs
% ---------------------------------------------------------------------
for varg = 1:length(varargin)

    if ischar(varargin{varg})
        % reserved keywords
        if strcmpi('unc',varargin{varg})
            inputargs = {varargin{varg-1}, 'uncorrected'};
            varargin{varg} = {}; varargin{varg - 1} = {};
        end
        if strcmpi('fdr',varargin{varg})
            inputargs = {varargin{varg-1}, 'fdr'};
            varargin{varg} = {}; varargin{varg - 1} = {};
        end

        if strcmpi('nodisplay',varargin{varg})
            do_display = 0;
            varargin{varg} = {};
        end

        if strcmpi('display',varargin{varg}) | strcmpi('display_results',varargin{varg})
            do_display = 0;
            varargin{varg} = {};
        end

        if strcmpi('robust',varargin{varg})
            do_robust = 1;
            varargin{varg} = {};
        end

        if strcmpi(varargin{varg}, 'AR')
            do_AR = true;
            % Check if the next argument is numeric (AR order)
            if (varg < length(varargin)) && isnumeric(varargin{varg+1})
                ar_order = varargin{varg+1};
                varargin{varg+1} = [];
            else
                ar_order = 1; % default AR(1) if none specified
            end
            varargin{varg} = [];
        end

        if strcmpi('brainony',varargin{varg}) | strcmpi('brain_is_predictor',varargin{varg}) %#ok<*OR2>
            brain_is_outcome = 0;
            varargin{varg} = {};
        end
        if strcmpi('nointercept',varargin{varg})
            do_intercept = 0;
            varargin{varg} = {};
        end
        if strcmpi('residual',varargin{varg})
            do_resid = 1;
            varargin{varg} = {};
        end

        if strcmpi('noverbose',varargin{varg})
            doverbose = false;
            varargin{varg} = {};
        end

        if strcmpi('names',varargin{varg}) | strcmpi('variable_names',varargin{varg})
            variable_names = varargin{varg + 1};
            varargin{varg} = {}; varargin{varg + 1} = {};
        end

        if strcmpi('C',varargin{varg}) | strcmpi('contrasts',varargin{varg})
            C = varargin{varg + 1};
            varargin{varg} = {}; varargin{varg + 1} = {};
        end

        if strcmpi('contrast_names',varargin{varg})
            contrast_names = varargin{varg + 1};
            varargin{varg} = {}; varargin{varg + 1} = {};
        end

        if strcmpi('analysis_name',varargin{varg})
            analysis_name = varargin{varg + 1};
            varargin{varg} = {}; varargin{varg + 1} = {};
        end

        if strcmpi('grandmeanscale',varargin{varg})
            grandmeanscale = true;
            varargin{varg} = {};
        end

        if strcmpi('covdat',varargin{varg})
            covdat = varargin{varg + 1};
            varargin{varg} = {}; varargin{varg + 1} = {};
        end

        if strcmpi('add_voxelwise_intercept', varargin{varg})
            add_voxelwise_intercept = true;
            varargin{varg} = {};
        end



    end % if ischar
end

if ~doverbose % add to pass into threshold( )
    inputargs{end+1} = 'noverbose';
end

if add_voxelwise_intercept & ~do_resid
    error('add_voxelwise_intercept = true but residuals not saved (do_resid = false). This is not compatible.')
end

% ---------------------------------------------------------------------
% Check Data and Diagnostics
% ---------------------------------------------------------------------

mywarnings = {};


% Check if fmri_data or image_vector
% ---------------------------------------------------------------------

if ~isa(obj,'fmri_data')
    error('obj input must be fmri_data object')
end

% Check data bit rate
% ---------------------------------------------------------------------

nuniquevals = length(unique(obj.dat(:)));

if nuniquevals < 2^10
    mywarnings{end+1} = sprintf('Number of unique values in dataset is low (%d, %3.2f bits), indicating possible restriction of bit rate. For comparison, Int16 has 65,536 unique values', nuniquevals, log2(nuniquevals));
end

% Check if Rank Deficient
% ---------------------------------------------------------------------

if rank(obj.X) < size(obj.X,2)
    mywarnings{end+1} = 'Warning:  obj.X is rank deficient.';
end

intercept_string = 'intercept is last';

% Check if Intercept is in model or specified for x_on_brain default
% ---------------------------------------------------------------------

if do_intercept && brain_is_outcome
    wh_int = intercept(obj.X, 'which');

    if isempty(wh_int)
        % add intercept and update wh_int (used later)
        if doverbose, mywarnings{end+1} = 'No intercept detected, adding intercept to last column of design matrix'; end
        X = intercept(obj.X, 'add');
        variable_names{end + 1} = 'Intercept';
        wh_int = intercept(X, 'which');
        obj.X(:, end+1) = repmat(1, size(obj.dat, 2), 1);
        X = obj.X;

    else
        intercept_string = sprintf('Intercept detected in column %1.0f of obj.X', wh_int);

        if doverbose, mywarnings{end+1} = intercept_string; end
        X = obj.X;

    end

else
    % No intercept or exogenous variables are outcome

    X = obj.X;
    intercept_string= 'No intercept included in model';

end

% Check of obj.X exists and is correct format
% ---------------------------------------------------------------------

if brain_is_outcome
    if isempty(obj.X)
        error('Make sure you include a design matrix in obj.X')
    end
    if size(obj.dat, 2) ~= size(obj.X, 1)
        error('obj.dat must have same number of columns as obj.X has rows.')
    end
    if isa(obj.X,'design_matrix')
        obj.X = obj.X.dat;
    end
else % Check if obj.Y exists and is in correct format if running brainony
    if isempty(obj.Y)
        error('Make sure you include a vector in obj.Y.')
    end
    if size(obj.dat, 2) ~= size(obj.Y, 1)
        error('obj.dat must have same number of columns as obj.Y has rows.')
    end
end

if do_resid && add_voxelwise_intercept

    wh_int = intercept(X, 'which'); % final value of which column is intercept
    if isempty(wh_int), error('no intercept in X. This is incompatible with add_voxelwise_intercept option.'), end

end

if do_resid && add_voxelwise_intercept

    wh_int = intercept(X, 'which'); % final value of which column is intercept
    if isempty(wh_int), error('no intercept in X. This is incompatible with add_voxelwise_intercept option.'), end

end

% Predictor centering
% ---------------------------------------------------------------------

m = mean(X);
wh_int = intercept(X, 'which');
m(wh_int) = [];
non_centered = abs(m) > 100 * eps;

% For effects-coded values [-1 1], ok, we want the intercept to reflect stats at average of 2 groups
iseffectcode = all(abs(X) == 1); % for each column
iseffectcode(wh_int) = [];

if any(non_centered & iseffectcode)

    mywarnings{end+1} = 'Warning:  Group sizes are unequal for effects-coded [1 -1] variable.';

end

if any(non_centered & ~iseffectcode)
    mywarnings{end+1} = 'Warning:  Predictors are not centered -- intercept is not interpretable as stats for average subject';
end

% Variance inflation
% ---------------------------------------------------------------------

vifs = getvif(X);

if any(vifs > 4)

    mywarnings{end+1} = 'Warning!!!  Design multicolinearity. Some regressors have variance inflation factors > 4. Check regression_results.diagnostics';

end

% Leverages
% ---------------------------------------------------------------------

H = X*pinv(X);
%H = X*inv(X'*X)*X'  will be identical if not rank deficient
leverages = diag(H);

if any(abs(zscore(leverages)) >= 3)
    mywarnings{end+1} = 'Warning!!!  Some observations have high leverage values relative to others, regression may be unstable. abs(z(leverage)) > 3';
end

% Names
% ---------------------------------------------------------------------

k = size(X, 2);
if length(variable_names) < k
    if ~isempty(variable_names), mywarnings{end+1} = 'Warning!!!  Too few variable names entered, less than size(X, 2). Names may be inaccurate.'; end % suppress warning if NO names entered

    for i = length(variable_names)+1:k %#ok<*FXUP>
        variable_names{i} = sprintf('R%d', i); %#ok<*AGROW>
    end
end

if length(variable_names) > k
    mywarnings{end+1} = 'Warning!!!  Too many variable names entered, more than size(X, 2). Names may be inaccurate.';

    variable_names = variable_names(1:k);
end

% Check contrasts
% ---------------------------------------------------------------------

if ~isempty(C) && ~(size(C, 1) == size(X, 2))
    % Do this *after* adding intercept to X if needed

    disp('Contrasts entered, but size(C, 1) does not equal size(X, 2).');
    if do_intercept
        disp('Must have a contrast entry for each column of X (including the intercept)');
    else
        disp('Must have a contrast entry for each column of X');
    end
    error('Quitting.')

end

% Contrast names
% ---------------------------------------------------------------------

if ~isempty(C)

    kc = size(C, 2);

    if length(contrast_names) < kc
        if ~isempty(contrast_names), mywarnings{end+1} = 'Warning!!!  Too few contrast names entered, less than size(C, 2). Names may be inaccurate.'; end % suppress warning if NO names entered

        for i = length(contrast_names)+1:kc
            contrast_names{i} = sprintf('Con%d', i); %#ok<*AGROW>
        end
    end

    if length(contrast_names) > kc
        mywarnings{end+1} = 'Warning!!!  Too many contrast names entered, more than size(C, 2). Names may be inaccurate.';

        contrast_names = contrast_names(1:kc);
    end

end

% Enforce proper shape
if ~iscolumn(variable_names), variable_names = variable_names'; end
if ~isrow(contrast_names), contrast_names = contrast_names'; end

% Enforce valid names: Eliminate special characters and leading numbers
[variable_names, namewarnings] = format_text_letters_only(variable_names, 'numbers', 'cleanup', 'squeeze', 'underscore_ok');
mywarnings = [mywarnings namewarnings];
if ~isempty(namewarnings), mywarnings{end+1} = 'Enter valid variable_names'; end

if ~isempty(C)

    [contrast_names, namewarnings] = format_text_letters_only(contrast_names, 'numbers', 'cleanup', 'squeeze', 'underscore_ok');
    mywarnings = [mywarnings namewarnings];
    if ~isempty(namewarnings), mywarnings{end+1} = 'Enter valid contrast_names'; end

end

% Data scaling and format
% ---------------------------------------------------------------------

% Enforce double-format (just in case)
obj.dat = double(obj.dat);

if grandmeanscale

    % Scale grand mean to 100; approximates what SPM and other packages do. Assumes mask and overall brain size are consistent across replicates (e.g., participants, in a first-level analysis
    % Scale run-wise, using obj.images_per_session to rescale each run to a grand mean of 100
    obj = obj.rescale('session_grand_mean_scaling_spm_style');

    % obj.dat = obj.dat .* 100 / nanmean(obj.dat(:));

end

% Check and set up for covdat (voxel-varying covariates)
% ---------------------------------------------------------------------
if ~isempty(covdat)

    % check robust enabled
    if ~do_robust
        warning('Covdat (voxel-varying covariates) requested, but robust flag not enabled. Covariates will be ignored (see Help)');
    end

    % if not in cell array, put it in cell array
    if iscell(covdat)
        % fine - nothing to do
    elseif isa(covdat, 'image_vector')
        covdat = {covdat};
    else
        error('Covdat must be a cell array of image_vectors or an image_vector')
    end

    for i=1:numel(covdat)
        % resample to obj space if needed
        if(compare_space(obj, covdat{i}))
            covdat{i} = resample_space(covdat{i}, obj);
        end

        % add covs to end of variable names
        if ~isempty(variable_names)
            variable_names{end+1} = sprintf('voxelwise_cov%d', i);
        end
    end

end


if doverbose

    fprintf('Analysis: %s\n', analysis_name);
    disp('----------------------------------');

    nowarnings = all(cellfun(@isempty, mywarnings));

    disp('Design matrix warnings:');
    disp('----------------------------------');
    if nowarnings
        disp('None')

    else
        disp(char(mywarnings{:}))

    end

    disp(' ');
end

% ---------------------------------------------------------------------
% Run Regression
% ---------------------------------------------------------------------

tic
s = warning;
warningwason = strcmp(s(1).state, 'on');
warning off

if brain_is_outcome
    % default is to regress obj.X on obj.dat (x on brain)
    % ---------------------------------------------------------------------

    % display info about regression
    if doverbose
        linestr = '______________________________________________________';
        disp(linestr);

        fprintf('Running regression: %3.0f voxels. Design: %3.0f obs, %3.0f regressors, %s\n', size(obj.dat, 1), size(X, 1), size(X, 2), intercept_string);
        if brain_is_outcome
            fprintf('\nPredicting exogenous variable(s) in obj.X using brain data as predictors, mass univariate');

        else % default
            fprintf('\nPredicting Brain Activity from obj.X, mass univariate');
        end
    end


    if do_AR
        if doverbose
            fprintf('\nRunning in GLS Mode with AR(%d)\n', ar_order);
        end

        v = size(obj.dat, 1);  % number of voxels
        [n, k] = size(X);

        % Pre-allocate outputs
        b = zeros(k, v);
        t = zeros(k, v);
        p = zeros(k, v);
        sigma = zeros(1, v);
        dfe = zeros(1, v);
        stderr = zeros(k, v);
        Phi_all = cell(1, v);

        % Pre-allocate contrast arrays if contrasts are provided
        if ~isempty(C)
            nc = size(C, 2);
            con_vals = zeros(nc, v);
            con_t = zeros(nc, v);
            con_p = zeros(nc, v);
        end

        % Loop over voxels and apply GLS via fit_gls2
        for i = 1:v
            y_voxel = obj.dat(i, :)';  % Voxel time-series as column vector
            [beta_i, t_i, pvals_i, con_vals_i, con_t_i, con_pvals_i, sigma_i, Phi_i, df_i, stderr_i, ~, ~] = ...
                fit_gls2(y_voxel, X, C, ar_order);

            b(:, i) = beta_i;
            t(:, i) = t_i;
            p(:, i) = pvals_i;
            sigma(:,i) = sigma_i;
            dfe(:,i) = df_i;
            stderr(:, i) = stderr_i;
            Phi_all{i} = Phi_i;

            if ~isempty(C)
                con_vals(:, i) = con_vals_i;
                con_t(:, i) = con_t_i;
                con_p(:, i) = con_pvals_i;
            end

            if doverbose && mod(i, 100) == 0
                fprintf('\b\b\b\b%03d%%', round(100 * i/v));
            end
        end % voxels

        % Compute residuals (if you need them for output)
        r = obj.dat' - X * b;

        if add_voxelwise_intercept

            r = r + b(wh_int, :);

        end

        % Continue with your existing robust or OLS branches

    elseif do_robust
        % non-AR
        % Robust - Regress stim/behavior(X) on brain (Y)
        %need to loop through voxels - Slow!

        if doverbose
            fprintf('\nRunning in Robust Mode ___%%');
        end

        v = size(obj.dat, 1);
        [n, k] = size(X);

        % multiple covariates
        if ~isempty(covdat)  % if we have covariate for each vox
            k = k + numel(covdat);
        end

        % Initialize outputs
        [b, t, stderr] = deal(zeros(k, v));
        p = ones(k, v);
        [dfe, sigma] = deal(zeros(1, v));


        for i = 1:v
            % For each voxel

            if ~isempty(covdat)  % voxel-varying covs: if we have different covariates for each vox

                Xi = X;
                for j=1:numel(covdat)
                    Xi(:,end+1) = covdat{j}.dat(i,:)';
                end

                [bb,stats] = robustfit(Xi, obj.dat(i,:)', 'bisquare', [], 'off');
                r(:,i) = stats.resid; % save voxelwise residual

                if add_voxelwise_intercept

                    r(:,i) = r(:,i) + bb(wh_int);

                end


            else % standard robust regression -- same covs for all vox

                [bb,stats] = robustfit(X, obj.dat(i,:)', 'bisquare', [], 'off');
            end

            b(:,i)=bb; %Betas
            t(:,i)=stats.t; %t-values
            p(:,i)=stats.p; %p-values
            dfe(:,i)=stats.dfe; %degrees of freedom
            stderr(:,i)=stats.se; %Standard Error
            sigma(:,i)=stats.robust_s; %robust estimate of sigma. LC not sure this is the best choice can switch to 'OLS_s','MAD_s', or 's'

            if doverbose && mod(i, 100) == 0
                fprintf('\b\b\b\b%03d%%', 100 * round(i/v))
            end

        end % voxels

        % if no vox varying covs, one-shot computation of residuals
        % residuals for vox varying covs are computed in loop above
        if isempty(covdat)
            r = obj.dat' - X*b; %residual

            if add_voxelwise_intercept

                r = r + b(wh_int, :);

            end

        end


    else % non-AR, non-robust
        % OLS - X predicting brain
        % - vectorized - Fast!

        if doverbose, fprintf('\nRunning in OLS Mode'); end

        % Estimate Betas in vector


        b = pinv(X) * obj.dat';

        % Error
        r = obj.dat' - X*b;

        % Residual variance
        [stderr, sigma] = get_std_errors(r, X);

        % Inference
        [t,dfe,p,sigma] = param_t_test(X,b,stderr,sigma);

        if add_voxelwise_intercept

            r = r + b(wh_int, :);

        end

    end % do-AR, do-robust, or do-OLS

else % brain is not outcome, it is predictor
    % Regress brain (X) on stim/behavior (Y) - loops through voxels, slow!
    % ---------------------------------------------------------------------

    if doverbose
        % display info about regression
        fprintf('regression > X: %3.0f voxels. Y: %3.0f obs, %3.0f regressors, %s\n', size(obj.dat, 1), size(obj.Y, 2), intercept_string);
        fprintf('\nPredicting obj.Y from Brain Activity');
    end

    if do_robust %need to loop through voxels - Slow!
        if doverbose
            fprintf('\nRunning in Robust Mode ___%%');
        end

        v = size(obj.dat, 1);
        n = size(obj.dat, 2);
        k = 2; % one predictor (brain) and an intercept

        % Initialize outputs
        [b, t, stderr] = deal(zeros(k, v));
        p = ones(k, v);
        [dfe, sigma] = deal(zeros(1, v));

        for i = 1:v
            % Create X from brain Data
            if do_intercept
                X = intercept(obj.dat(i,:)','add');
            else
                X = obj.dat(i,:)';
            end

            [bb,stats] = robustfit(X, obj.Y, 'bisquare', [], 'off');

            b(:,i)=bb; %Betas
            t(:,i)=stats.t; %t-values
            p(:,i)=stats.p; %p-values
            dfe(:,i)=stats.dfe; %degrees of freedom
            stderr(:,i)=stats.se; %Standard Error
            sigma(:,i)=stats.robust_s; %robust estimate of sigma. LC not sure this is the best choice can switch to 'OLS_s','MAD_s', or 's'
            r(:,i) = obj.Y - X * b(:,i); %residual

            if doverbose && mod(i, 100) == 0
                fprintf('\b\b\b\b%03d%%', 100 * round(i/v))
            end

        end % voxel

    else % non-robust
        % ---------------------------------------------------------------------
        %OLS -- Regress brain on Y

        if doverbose, fprintf('\nRunning in OLS Mode'); end

        for i = 1:size(obj.dat,1)

            % Create X from brain Data
            if do_intercept
                X = intercept(obj.dat(i,:)','add');
            else
                X = obj.dat(i,:)';
            end

            % Estimate Betas in vector
            b(:,i) = pinv(X) * obj.Y;

            % Error
            r(:,i) = obj.Y - X * b(:,i);

            %             sigma(i) = std(r(:,i)); % wrong
            %             stderr(:,i) = ( diag(inv(X' * X)) .^ .5 ) * sigma(i);  % params x voxels matrix of std. errors
        end

        % Residual variance - can do this at end because X is always the same size
        [stderr, sigma] = get_std_errors(r, X);

        % Inference
        [t,dfe,p,sigma] = param_t_test(X,b,stderr,sigma);

    end % dorobust

end % brain_is_outcome or not

stop = toc;
if doverbose, fprintf('\nModel run in %d minutes and %.2f seconds\n',floor(stop/60),rem(stop,60)); end


% ---------------------------------------------------------------------
% Contrasts
% ---------------------------------------------------------------------

if ~isempty(C)

    con_vals = C' * b;

    % Contrast STE
    con_ste = diag(C' * inv(X' * X) * C) .^ .5 * sigma;

    [con_t, ~, con_p] = param_t_test(X, con_vals, con_ste, sigma);

end


% ---------------------------------------------------------------------
% Create Output - out = regression_results
% ---------------------------------------------------------------------

out = struct;

out.analysis_name = analysis_name;

out.input_parameters = struct( ...
    'brain_is_predictor', brain_is_outcome, 'do_robust', do_robust, 'grandmeanscale', grandmeanscale, ...
    'do_intercept', do_intercept, ...
    'do_resid', do_resid, 'doverbose', doverbose, 'do_display', do_display, 'covdat', covdat);

out.input_parameters.initial_statistical_threshold = inputargs;

out.input_image_metadata.source_notes = obj.source_notes;
out.input_image_metadata.history = obj.history;
out.input_image_metadata.image_names = obj.image_names;
out.input_image_metadata.fullpath = obj.fullpath;

% design, contrasts, and diagnostics

out.X = X;
out.variable_names = variable_names;
out.C = C;
out.contrast_names = contrast_names;

out.contrast_summary_table = table();

if ~isempty(C)
    % Contrast summary table

    for i = 1:size(C, 2)
        out.contrast_summary_table(:, i) = table(C(:, i));
    end

    out.contrast_summary_table.Properties.RowNames = variable_names;
    out.contrast_summary_table.Properties.VariableNames = contrast_names;

    if doverbose
        fprintf('\nSummary of conditions and contrasts\n%s\n', linestr);
        disp(out.contrast_summary_table);
    end

end

out.diagnostics = struct('Variance_inflation_factors', vifs, 'Leverages', leverages);
out.warnings = mywarnings;

% Create objects
if doverbose
    fprintf('\nCreating beta and t images, thresholding t images\n%s\n', linestr);
end

% Betas
out.b = statistic_image;
out.b.type = 'Beta';
out.b.p = p';
out.b.ste = stderr';
out.b.N = n;
out.b.dat = b';
out.b.dat_descrip = sprintf('Beta Values from regression, intercept is column %d', wh_int);
out.b.volInfo = obj.volInfo;
out.b.removed_voxels = obj.removed_voxels;
out.b.removed_images = false;  % this image does not have the same dims as the original dataset
out.b.image_labels = variable_names;

out.b = enforce_variable_types(out.b);

% - beta and contrast images are unthresholded, t images from both are thresholded
% if doverbose, fprintf('Thresholding b images at %3.6f %s\n', inputargs{1}, inputargs{2}); end
% out.b = threshold(out.b, inputargs{:}, 'noverbose'); % Threshold image

% T stats
out.t = statistic_image;
out.t.type = 'T';
out.t.p = p';
out.t.ste = stderr';
out.t.N = n;
out.t.dat = t';
out.t.dat_descrip = sprintf('t-values from regression, intercept is column %d', wh_int);
out.t.volInfo = obj.volInfo;
out.t.removed_voxels = obj.removed_voxels;
out.t.removed_images = false;  % this image does not have the same dims as the original dataset
out.t.image_labels = variable_names;

% out.t = enforce_variable_types(out.t); % The remove_empty in this can cause some problems for regress

if doverbose
    fprintf('Thresholding t images at %3.6f %s\n', inputargs{1}, inputargs{2});
end
out.t = threshold(out.t, inputargs{:}); %Threshold image

% DF as fmri_data
out.df = obj;
out.df.dat = dfe';
out.df.dat_descrip = sprintf('Degrees of Freedom');

% Sigma as fmri_data
out.sigma = obj;
out.sigma.dat = sigma';
out.sigma.dat_descrip = sprintf('Sigma from Regression');

% Residual as fmri_data
if do_resid
    out.resid = obj;
    out.resid.dat = r';
    out.resid.dat_descrip = sprintf('Residual from Regression');
    out.resid.history{end + 1} = 'Removed covariates via regression';
end

if ~isempty(C)

    if doverbose, fprintf('\nCreating contrast and t images and thresholding t images\n%s\n', linestr); end

    % Contrast values
    out.contrast_images = statistic_image;
    out.contrast_images.type = 'Contrast';
    out.contrast_images.p = con_p';
    out.contrast_images.ste = stderr';
    out.contrast_images.N = n;
    out.contrast_images.dat = con_vals';
    out.contrast_images.dat_descrip = 'Contrast Values from regression';
    out.contrast_images.volInfo = obj.volInfo;
    out.contrast_images.removed_voxels = obj.removed_voxels;
    out.contrast_images.removed_images = false;  % this image does not have the same dims as the original dataset
    out.contrast_images.image_labels = contrast_names;

    out.contrast_images = enforce_variable_types(out.contrast_images);

    %     if doverbose, fprintf('Thresholding contrast images at %3.6f %s\n', inputargs{1}, inputargs{2}); end
    %     out.contrast_images = threshold(out.contrast_images, inputargs{:}, 'noverbose'); % Threshold image

    % T stats
    out.con_t = statistic_image;
    out.con_t.type = 'T';
    out.con_t.p = con_p';
    out.con_t.ste = con_ste';
    out.con_t.N = n;
    out.con_t.dat = con_t';
    out.con_t.dat_descrip = sprintf('t-values from regression, intercept is column %d', wh_int);
    out.con_t.volInfo = obj.volInfo;
    out.con_t.removed_voxels = obj.removed_voxels;
    out.con_t.removed_images = false;  % this image does not have the same dims as the original dataset
    out.con_t.image_labels = contrast_names;

    out.con_t = enforce_variable_types(out.con_t);

    if doverbose
        fprintf('Thresholding t images at %3.6f %s\n', inputargs{1}, inputargs{2});

    end
    out.con_t = threshold(out.con_t, inputargs{:}); %Threshold image

end


% ---------------------------------------------------------------------
% Plot Results
% ---------------------------------------------------------------------
if k < 10 && do_display

    orthviews(out.t);

elseif do_display && doverbose
    disp('Warning: No display because >= 10 images.');

end

% reset warning state if needed
if warningwason
    warning on
end

% ---------------------------------------------------------------------
% Subfunctions
% ---------------------------------------------------------------------



    function [stderr, sigma] = get_std_errors(r, X)

        [n, k] = size(X);
        v = size(obj.dat, 1);

        % We want diag(r' * r), but matrix size can be large
        % std(r) does not account for k parameters used, so is incorrect
        % residual std. not sqrt(var(resid)) -- we must account for k params used
        for i = 1:v
            sigma(1, i) = (r(:, i)' * r(:, i) ./ (n - k)) .^ .5;
        end

        stderr = ( diag(inv(X' * X)) .^ .5 ) * sigma;  % params x voxels matrix of std. errors

    end


    function [t,dfe,p,sigma] = param_t_test(X,b,stderr,sigma)
        % test whether parameter is significantly different from zero
        %
        % Inputs:
        % X:        design matrix
        % b:        beta values
        % stderr:   standard error of beta estimate
        % sigma:    standard deviation of residual
        %
        % Returns:
        % t:        t-value
        % dfe:      degrees of freedom
        % p:        p-value
        % sigma:    standard deviation of residual

        [n, k] = size(X);
        t = b ./ stderr;
        dfe = n - k;
        p = 2 * (1 - tcdf(abs(t), dfe));

        sigma(sigma == 0) = Inf;
        t(isnan(t)) = 0;
        p(isnan(p)) = 0;
        dfe = repmat(dfe,1,size(t,2));
    end

end % Main Function
