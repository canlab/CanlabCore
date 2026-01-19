function out = fitlme_voxelwise(obj, tbl, spec, varargin)
% Voxelwise linear mixed-effects (LME) model using fitlme.
%
% - Uses obj.dat as brain data: [voxels x observations]
% - Uses tbl as design/covariate table: one row per observation
% - Third argument (spec) can be:
%       (a) A fitlme formula string, e.g.:
%           'Y ~ 1 + Run + Condition + (1 + Run | Subject)'
%       (b) A role table with variables:
%           Name : variable name in tbl
%           Role : 'group' | 'fixed' | 'mixed' | 'random_only'
%           In that case, the function builds the formula automatically.
%
% - Runs the same mixed model at each voxel
% - Returns beta/t maps as statistic_image objects
% - Supports fixed-effect contrasts (matrix C) like regress.m
%
% :Usage (formula mode):
% ::
%
%   out = fitlme_voxelwise(obj, tbl, ...
%       'Y ~ 1 + Run + Condition + (1 + Run | Subject)');
%
%   out = fitlme_voxelwise(obj, tbl, ...
%       'Y ~ 1 + Run + Condition + (1 + Run | Subject)', ...
%       .005, 'unc', 'C', C, 'contrast_names', {...});
%
% :Usage (role-table mode):
% ::
%
%   % tbl has columns: Subject, Baseline, Run, C1, C2
%   role_tbl = table( ...
%       {'Subject'; 'Baseline'; 'Run'; 'C1'; 'C2'}, ...
%       {'group';   'fixed';    'mixed'; 'mixed'; 'mixed'}, ...
%       'VariableNames', {'Name','Role'});
%
%   out = fitlme_voxelwise(obj, tbl, role_tbl, ...
%       .001, 'unc', 'analysis_name', 'LME_with_roles');
%
%   % This internally builds:
%   %   Y ~ 1 + Baseline + Run + C1 + C2 + (1 + Run + C1 + C2 | Subject)
%
% :Inputs:
%
%   **obj:**
%       fmri_data (recommended) or image_vector with:
%       obj.dat  = [nVox x nObs]
%
%   **tbl:**
%       MATLAB table with nObs rows.
%       Does NOT need a 'Y' column; that is added internally per voxel.
%
%   **spec:**
%       EITHER:
%         - fitlme formula string ('Y ~ ...')
%         - OR role table with Name/Role columns
%
% :Optional Inputs (similar spirit to regress.m):
%
%   **pthr, 'unc' or 'fdr':**
%       Threshold for t-maps / contrast t-maps.
%       Example: .005, 'unc'  or  .05, 'fdr'
%       If omitted, default is p = .001, 'uncorrected'.
%
%   **'C', C:**
%       Contrast matrix over *fixed effects*.
%       Size: [nFixed x nContrasts], where nFixed is number of fixed effects
%       (including intercept) returned by fitlme.
%
%   **'contrast_names', {cellstr}:**
%       Names for each contrast (columns of C).
%
%   **'analysis_name', string:**
%       Name/description stored in out.analysis_name.
%
%   **'doparallel', logical:**
%       Use parfor over voxels (if Parallel Toolbox available).
%       Default = false.
%
%   **'noverbose':**
%       Suppress verbose text output.
%
%   **'residual':**
%       Save residual time series per voxel as out.resid (fmri_data).
%
%   **'fitmethod', 'REML' or 'ML':**
%       Passed to fitlme (default 'REML').
%
%   Unsupported options ('robust', 'AR', 'brainony', 'covdat', etc.)
%   will be ignored with a warning.
%
% :Outputs (CANlab-style):
%
%   **out.b:**
%       statistic_image of fixed-effect betas
%         - out.b.dat         = [nVox x nFixed]
%         - out.b.image_labels = fixed effect names
%
%   **out.t:**
%       statistic_image of fixed-effect t-values (thresholded)
%
%   **out.df:**
%       fmri_data: df per voxel (scalar DFE from fitlme, first effect)
%
%   **out.sigma:**
%       fmri_data: residual std dev per voxel
%
%   **out.contrast_images (if C provided):**
%       statistic_image of contrast beta values
%
%   **out.con_t (if C provided):**
%       statistic_image of contrast t-values (thresholded)
%
%   **out.resid (optional):**
%       fmri_data of residuals [nVox x nObs] (transposed at end)
%
%   **out.fixed_names:**
%       cellstr of fixed-effect names from fitlme
%
%   **out.C, out.contrast_names:**
%       contrast matrix and names (if given)
%
%   **out.formula, out.input_parameters, out.warnings:**
%       metadata about the analysis
%
% Examples:
% -----------------------------------------------
% % Run on BMRK3 pain dataset, 6 temperatures for each of 33 subjects
% % Random intercept/random slope, random effect of temperature
%
% fmri_data_file = which('bmrk3_6levels_pain_dataset.mat');
% load(fmri_data_file)
%
% % Create metadata_table with images x variables matrix with names
% t = table(image_obj.additional_info.subject_id, image_obj.additional_info.temperatures, 'VariableNames', {'subject_id', 'temperature'});
% image_obj.metadata_table = t;
%
% % Run fitlme on every voxel using Wilkinson notation to refer to the table:
% out = fitlme_voxelwise(image_obj, image_obj.metadata_table, 'Y ~ 1 + temperature + (1 + temperature | subject_id)', .005, 'unc');
% t_obj_temp_fitlme = get_wh_image(out.t, 2);
% 
% % Compare t-values for fitlme against two-stage summary statistics approach
% % First, create an object with a linear contrast across temps for each person:
% C = [-5 -3 -1 1 3 5]';
% C_subj = repmat({C}, 1, 33);
% C_subj_matrix = blkdiag(C_subj{:});
% contrast_dat = image_obj.dat * C_subj_matrix;
% contrast_obj = image_obj;
% contrast_obj.dat = contrast_dat;
% contrast_obj.removed_images = false(33, 1);
% contrast_obj.dat_descrip = '.dat contains one linear contrast across temperatures for each of 33 subjects';
%
% % Fit 2SSS approach and get t (one sample t-test)
% t_obj_temp_2sss = ttest(contrast_obj);
% t_obj_temp_2sss = threshold(t_obj_temp_2sss, 0.005, 'unc');
% create_figure; axis off; montage(t_obj_temp_2sss, 'onerow')
%
% % Compare t-maps in a scatterplot
% figure; plot(t_obj_temp_2sss.dat, t_obj_temp_fitlme.dat, '.');
% hold on; plot([-10 10], [-10 10], '--', 'Color', 'k');
% ylabel('t-values after normalization'); xlabel('t-values before normalization');
%
% or:
% h = image_scatterplot(t_obj_temp_2sss, t_obj_temp_fitlme, 'pvaluebox', 0.001, 'colorpoints');

% 
% -------------------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------------------

if nargin < 3 || isempty(spec)
    error('fitlme_voxelwise:TooFewInputs', ...
        'Usage: out = fitlme_voxelwise(obj, tbl, formula_or_role_tbl, [...])');
end

if ~istable(tbl)
    error('tbl must be a MATLAB table.');
end

if ~isa(obj, 'fmri_data') && ~isa(obj, 'image_vector')
    error('obj must be fmri_data or image_vector (with .dat and .volInfo).');
end

% obj is an object, so use isprop not isfield:
if ~isprop(obj, 'dat') || isempty(obj.dat)
    error('obj.dat is missing or empty.');
end

Ymat = double(obj.dat);            % [nVox x nObs], use double
[nVox, nObs] = size(Ymat);

if height(tbl) ~= nObs
    error('Number of observations in obj.dat (columns) must match height(tbl).');
end

% -------------------------------------------------------------------------
% Build formula from spec (string OR role-table)
% -------------------------------------------------------------------------

if ischar(spec) || isstring(spec)
    formula = char(spec);  % user supplied formula
elseif istable(spec)
    role_tbl = spec;
    
    % Validate role_tbl
    if ~all(ismember({'Name','Role'}, role_tbl.Properties.VariableNames))
        error('Role table must have variables: Name, Role.');
    end
    
    names = cellstr(role_tbl.Name);
    roles = lower(cellstr(role_tbl.Role));
    
    % Check names exist in tbl
    if ~all(ismember(names, tbl.Properties.VariableNames))
        missing = setdiff(names, tbl.Properties.VariableNames);
        error('Variables in role_tbl not found in tbl: %s', strjoin(missing, ', '));
    end
    
    groupvar              = names(strcmp(roles, 'group'));
    fixed_only_vars       = names(strcmp(roles, 'fixed'));
    random_and_fixed_vars = names(strcmp(roles, 'mixed'));
    random_only_vars      = names(strcmp(roles, 'random_only'));
    
    if numel(groupvar) ~= 1
        error('Exactly one variable must have Role == ''group'' (grouping factor).');
    end
    groupvar = groupvar{1};
    
    % Fixed part: 1 + fixed_only + random_and_fixed
    fixed_terms = {'1'};
    if ~isempty(fixed_only_vars)
        fixed_terms = [fixed_terms, fixed_only_vars(:)'];
    end
    if ~isempty(random_and_fixed_vars)
        fixed_terms = [fixed_terms, random_and_fixed_vars(:)'];
    end
    fixed_terms = unique(fixed_terms, 'stable');
    fixed_str = strjoin(fixed_terms, ' + ');
    
    % Random part: (1 + random_slopes | groupvar)
    random_slopes = unique([random_and_fixed_vars(:); random_only_vars(:)], 'stable');
    random_terms  = [{'1'}; random_slopes];
    random_str    = sprintf('(%s | %s)', strjoin(random_terms, ' + '), groupvar);
    
    formula = sprintf('Y ~ %s + %s', fixed_str, random_str);
else
    error('Third input must be a formula string or a Name/Role table.');
end

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------

pthr         = 0.001;
thresh_type  = 'uncorrected';  % pass into threshold; compatible with 'unc' / 'fdr'
do_display   = false;
doverbose    = true;
doparallel   = false;
do_resid     = false;
fitmethod    = 'REML';
analysis_name = '';

C = [];
contrast_names = {};

mywarnings = {};

% -------------------------------------------------------------------------
% Parse varargin (similar flavor to regress.m)
% -------------------------------------------------------------------------

i = 1;
while i <= length(varargin)
    arg = varargin{i};
    
    if isnumeric(arg) && i < length(varargin) && ischar(varargin{i+1})
        % threshold pair: pthr, 'unc' or 'fdr'
        pthr = arg;
        tstr = lower(varargin{i+1});
        if strcmp(tstr, 'unc')
            thresh_type = 'uncorrected';
        elseif strcmp(tstr, 'fdr')
            thresh_type = 'fdr';
        else
            thresh_type = tstr;  % allow passing through
        end
        i = i + 2;
        continue;
    end
    
    if ischar(arg)
        switch lower(arg)
            
            case {'display','display_results'}
                do_display = true;
                i = i + 1;
                
            case 'nodisplay'
                do_display = false;
                i = i + 1;
                
            case 'noverbose'
                doverbose = false;
                i = i + 1;
                
            case 'doparallel'
                doparallel = true;
                i = i + 1;
                
            case 'residual'
                do_resid = true;
                i = i + 1;
                
            case 'fitmethod'
                if i == length(varargin)
                    error('Missing value after ''fitmethod''.');
                end
                fitmethod = varargin{i+1};
                i = i + 2;
                
            case {'analysis_name'}
                if i == length(varargin)
                    error('Missing value after ''analysis_name''.');
                end
                analysis_name = varargin{i+1};
                i = i + 2;
                
            case {'contrast_names'}
                if i == length(varargin)
                    error('Missing value after ''contrast_names''.');
                end
                contrast_names = varargin{i+1};
                i = i + 2;
                
            case {'c','contrasts'}
                if i == length(varargin)
                    error('Missing value after ''C'' / ''contrasts''.');
                end
                C = varargin{i+1};
                i = i + 2;
                
            % --- Unsupported regress options: warn & ignore ---
            case {'robust','brainony','brain_is_predictor','grandmeanscale','covdat','ar'}
                mywarnings{end+1} = sprintf( ...
                    'Option ''%s'' is not implemented for fitlme_voxelwise and will be ignored.', arg);
                i = i + 1;
                
            otherwise
                warning('Unknown option "%s" ignored.', arg);
                i = i + 1;
        end
        
    else
        warning('Unrecognized input at position %d; ignoring.', i);
        i = i + 1;
    end
end

if doverbose
    fprintf('fitlme_voxelwise: pthr = %.4g, type = %s\n', pthr, thresh_type);
    fprintf('  doparallel = %d, fitmethod = %s\n', doparallel, fitmethod);
    fprintf('  formula: %s\n', formula);
end

% -------------------------------------------------------------------------
% Initial model: get fixed-effect structure from first voxel
% -------------------------------------------------------------------------

tbl_work = tbl;
tbl_work.Y = Ymat(1, :)';

try
    lme0 = fitlme(tbl_work, formula, 'FitMethod', fitmethod);
catch ME
    error('Initial fitlme failed at first voxel: %s', ME.message);
end

[beta0, ~, stats0] = fixedEffects(lme0, 'DFMethod', 'Residual'); %#ok<ASGLU>
fixed_names = lme0.CoefficientNames(:)';  % includes intercept if present
nFixed      = numel(beta0);

if doverbose
    fprintf('Fixed effects (%d):\n', nFixed);
    disp(fixed_names);
end

% -------------------------------------------------------------------------
% Validate contrasts C
% -------------------------------------------------------------------------

if ~isempty(C)
    [kC, nC] = size(C);
    if kC ~= nFixed
        error('Contrast matrix C must have size [nFixed x nContrasts], here %d x nC, but nFixed = %d.', kC, nFixed);
    end
    
    if ~isempty(contrast_names)
        if numel(contrast_names) ~= nC
            warning('Length of contrast_names does not match number of contrasts; trimming or padding.');
            contrast_names = contrast_names(:)';
            if numel(contrast_names) < nC
                for ii = numel(contrast_names)+1:nC
                    contrast_names{ii} = sprintf('Con%d', ii);
                end
            else
                contrast_names = contrast_names(1:nC);
            end
        end
    else
        contrast_names = arrayfun(@(ii) sprintf('Con%d', ii), 1:nC, 'UniformOutput', false);
    end
else
    nC = 0;
end

% -------------------------------------------------------------------------
% Preallocate voxelwise containers
% -------------------------------------------------------------------------

betas = NaN(nFixed, nVox);
tvals = NaN(nFixed, nVox);
pvals = NaN(nFixed, nVox);
SE    = NaN(nFixed, nVox);
DFmat = NaN(nFixed, nVox);
sigma = NaN(1, nVox);
mask  = false(nVox, 1);

if do_resid
    R = NaN(nObs, nVox);  % residual matrix [obs x vox]
end

if nC > 0
    con_vals = NaN(nC, nVox);
    con_t    = NaN(nC, nVox);
    con_p    = NaN(nC, nVox);
    con_se   = NaN(nC, nVox);
    con_df   = NaN(nC, nVox);
else
    % keep 0 x nVox so indexing works even if never used
    con_vals = NaN(0, nVox);
    con_t    = NaN(0, nVox);
    con_p    = NaN(0, nVox);
    con_se   = NaN(0, nVox);
    con_df   = NaN(0, nVox);
end

% -------------------------------------------------------------------------
% Main voxelwise loop
% -------------------------------------------------------------------------

if doverbose
    fprintf('Running voxelwise LME: %d voxels, %d observations.\n', nVox, nObs);
    tic
end

if doparallel
    if doverbose, fprintf('Using PARFOR over voxels.\n'); end
    
    parfor v = 1:nVox
        [betas(:,v), tvals(:,v), pvals(:,v), SE(:,v), DFmat(:,v), ...
         sigma(1,v), res_v, ...
         con_vals(:,v), con_t(:,v), con_p(:,v), con_se(:,v), con_df(:,v), mask(v)] ...
         = fit_lme_one_voxel( ...
                        Ymat(v,:)', tbl, formula, fitmethod, C, nC, nFixed);
        if do_resid
            R(:,v) = res_v;
        end
    end
else
    if doverbose
            fprintf('  Serial run of %3.0f voxels: %3.0f%%', nVox, 0)
    end
            
    for v = 1:nVox
        
        if doverbose && rem(v, round(nVox/100)) == 0
            fprintf('\b\b\b\b%3.0f%%', round(100*v/nVox));
        end
        
        [betas(:,v), tvals(:,v), pvals(:,v), SE(:,v), DFmat(:,v), ...
         sigma(1,v), res_v, ...
         con_vals(:,v), con_t(:,v), con_p(:,v), con_se(:,v), con_df(:,v), ...
         mask(v)] = fit_lme_one_voxel( ...
                        Ymat(v,:)', tbl, formula, fitmethod, C, nC, nFixed);
        if do_resid
            R(:,v) = res_v;
        end
    end
end

if doverbose
    fprintf('Voxelwise LME finished. Successful fits: %d / %d voxels.\n', sum(mask), nVox);
    toc
end

% -------------------------------------------------------------------------
% Build CANlab-style outputs
% -------------------------------------------------------------------------

out = struct;
out.analysis_name    = analysis_name;
out.formula          = formula;
out.fixed_names      = fixed_names;
out.C                = C;
out.contrast_names   = contrast_names;
out.mask             = mask;
out.input_parameters = struct( ...
    'pthr', pthr, 'thresh_type', thresh_type, ...
    'doparallel', doparallel, 'fitmethod', fitmethod);
out.warnings         = mywarnings;

% --- df image (take first fixed-effect df as representative per voxel) ---
df_scalar        = DFmat(1, :);   % [1 x nVox]
df_img           = obj;
df_img.dat       = df_scalar';
df_img.dat_descrip = 'Degrees of freedom (first fixed effect DFE)';
out.df           = df_img;

% --- sigma image ---
sig_img          = obj;
sig_img.dat      = sigma';
sig_img.dat_descrip = 'Residual std dev from LME per voxel';
out.sigma        = sig_img;

% --- Beta maps as statistic_image ---
b_img                    = statistic_image;
b_img.type               = 'Beta (LME)';
b_img.dat                = betas';                 % [nVox x nFixed]
b_img.p                  = pvals';                % [nVox x nFixed]
b_img.ste                = SE';                   % [nVox x nFixed]
b_img.N                  = nObs;
b_img.volInfo            = obj.volInfo;
b_img.removed_voxels     = obj.removed_voxels;
b_img.removed_images     = false;
b_img.image_labels       = fixed_names;
b_img.dat_descrip        = 'Fixed-effect betas from voxelwise LME';
b_img                    = enforce_variable_types(b_img);
out.b                    = b_img;

% --- T maps as statistic_image ---
t_img                    = statistic_image;
t_img.type               = 'T (LME fixed effects)';
t_img.dat                = tvals';     % [nVox x nFixed]
t_img.p                  = pvals';     % [nVox x nFixed]
t_img.ste                = SE';        % [nVox x nFixed]
t_img.N                  = nObs;
t_img.volInfo            = obj.volInfo;
t_img.removed_voxels     = obj.removed_voxels;
t_img.removed_images     = false;
t_img.image_labels       = fixed_names;
t_img.dat_descrip        = 'Fixed-effect t-values from voxelwise LME';
t_img                    = enforce_variable_types(t_img);

if doverbose
    fprintf('Thresholding fixed-effect t-maps at p = %.4g (%s)\n', pthr, thresh_type);
end
t_img = threshold(t_img, pthr, thresh_type);
out.t = t_img;

% --- Residuals as fmri_data (optional) ---
if do_resid
    resid_img             = obj;
    resid_img.dat         = R';   % [nVox x nObs]
    resid_img.dat_descrip = 'Residuals from voxelwise LME';
    resid_img.history{end+1} = 'Residuals from fitlme_voxelwise';
    out.resid             = resid_img;
end

% --- Contrasts (if any) ---
if nC > 0
    % Contrast beta
    con_img                    = statistic_image;
    con_img.type               = 'Contrast (LME fixed effects)';
    con_img.dat                = con_vals';          % [nVox x nC]
    con_img.p                  = con_p';             % [nVox x nC]
    con_img.ste                = con_se';            % [nVox x nC]
    con_img.N                  = nObs;
    con_img.volInfo            = obj.volInfo;
    con_img.removed_voxels     = obj.removed_voxels;
    con_img.removed_images     = false;
    con_img.image_labels       = contrast_names;
    con_img.dat_descrip        = 'Contrast beta values from voxelwise LME';
    con_img                    = enforce_variable_types(con_img);
    out.contrast_images        = con_img;
    
    % Contrast t
    con_t_img                  = statistic_image;
    con_t_img.type             = 'T (LME contrasts)';
    con_t_img.dat              = con_t';             % [nVox x nC]
    con_t_img.p                = con_p';             % [nVox x nC]
    con_t_img.ste              = con_se';            % [nVox x nC]
    con_t_img.N                = nObs;
    con_t_img.volInfo          = obj.volInfo;
    con_t_img.removed_voxels   = obj.removed_voxels;
    con_t_img.removed_images   = false;
    con_t_img.image_labels     = contrast_names;
    con_t_img.dat_descrip      = 'Contrast t-values from voxelwise LME';
    con_t_img                  = enforce_variable_types(con_t_img);
    
    if doverbose
        fprintf('Thresholding contrast t-maps at p = %.4g (%s)\n', pthr, thresh_type);
    end
    con_t_img = threshold(con_t_img, pthr, thresh_type);
    out.con_t = con_t_img;
end

% -------------------------------------------------------------------------
% Display (optional)
% -------------------------------------------------------------------------
if do_display
    try
        orthviews(out.t);
    catch
        warning('orthviews failed; skipping display.');
    end
end

end % main function


% =====================================================================
% Helper: fit one voxel's LME and (optionally) contrasts
% =====================================================================
function [b, t, p, se, df, sigma, res_v, ...
          con_vals, con_t, con_p, con_se, con_df, ok] = ...
          fit_lme_one_voxel(Yvec, tbl, formula, fitmethod, C, nC, nFixed)

% Initialize with correct sizes so assignment in parfor works
b       = NaN(nFixed, 1);
t       = NaN(nFixed, 1);
p       = NaN(nFixed, 1);
se      = NaN(nFixed, 1);
df      = NaN(nFixed, 1);
sigma   = NaN;
res_v   = NaN(height(tbl), 1);
ok      = false;

if nC > 0
    con_vals = NaN(nC,1);
    con_t    = NaN(nC,1);
    con_p    = NaN(nC,1);
    con_se   = NaN(nC,1);
    con_df   = NaN(nC,1);
else
    con_vals = NaN(0,1);
    con_t    = NaN(0,1);
    con_p    = NaN(0,1);
    con_se   = NaN(0,1);
    con_df   = NaN(0,1);
end

% All-NaN Y -> skip (returns NaNs with correct sizes)
if all(isnan(Yvec))
    return;
end

tbl_local      = tbl;
tbl_local.Y    = Yvec;

try
    lme = fitlme(tbl_local, formula, 'FitMethod', fitmethod);
    
    [beta_hat, CovB, stats] = fixedEffects(lme, 'DFMethod', 'Residual');
    
    % overwrite initialized NaNs with real values
    k = numel(beta_hat);
    b(1:k)   = beta_hat;
    t(1:k)   = stats.tStat(:);
    p(1:k)   = stats.pValue(:);
    se(1:k)  = stats.SE(:);
    
    % DFE may be scalar or per-parameter; handle both
    if isscalar(stats.DFE)
        df(:) = stats.DFE;
    else
        df(1:k) = stats.DFE(:);
    end
    
    res_v = residuals(lme, 'ResidualType', 'Raw');
    sigma = std(res_v);
    
    % Contrasts (over fixed effects)
    if nC > 0
        for ci = 1:nC
            Li = C(:,ci);                 % [nFixed x 1]
            cval = Li' * beta_hat;        % scalar
            cse  = sqrt( Li' * CovB * Li );
            if cse == 0 || isnan(cse)
                ct = NaN;
                cp = NaN;
            else
                ct = cval / cse;
                dfe_con = stats.DFE;      % usually scalar
                if isscalar(dfe_con)
                    cp = 2 * (1 - tcdf(abs(ct), dfe_con));
                else
                    dfe_con = dfe_con(1);
                    cp = 2 * (1 - tcdf(abs(ct), dfe_con));
                end
            end
            con_vals(ci,1) = cval;
            con_se(ci,1)   = cse;
            con_t(ci,1)    = ct;
            con_p(ci,1)    = cp;
            if isscalar(stats.DFE)
                con_df(ci,1) = stats.DFE;
            else
                con_df(ci,1) = stats.DFE(1);
            end
        end
    end
    
    ok = true;
catch
    ok = false;
end

end
