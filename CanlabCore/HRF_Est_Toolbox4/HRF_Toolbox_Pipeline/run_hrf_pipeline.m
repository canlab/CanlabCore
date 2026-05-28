function results = run_hrf_pipeline(fmri_nii, events_tsv, varargin)
%RUN_HRF_PIPELINE End-to-end HRF estimation from 4D fMRI + BIDS events.
%   results = RUN_HRF_PIPELINE(fmri_nii, events_tsv, ...)
%
% Required inputs
%   fmri_nii    - path to 4D fMRI NIfTI file.
%   events_tsv  - path to BIDS-style events.tsv file with columns:
%                 onset, duration, trial_type.
%
% Name-value options
%   'TR'              : repetition time in seconds (default: from NIfTI header)
%   'MaskNii'         : optional mask NIfTI; if omitted uses whole-brain mean
%   'Conditions'      : cellstr of trial_type names to fit (default: all)
%   'WindowSeconds'   : HRF estimation window, seconds (default 30)
%   'Models'          : subset of {'logit','fir','sfir','canonical','spline','nlgamma'}
%                       (default all)
%   'WholeBrainMode'  : 'auto', a model name, or cellstr of model names.
%                       'auto' writes one whole-brain output per requested
%                       linear model: fir, sfir, canonical, and spline.
%   'OutputMat'       : optional .mat path to save results
%
% Output
%   results struct with fields:
%     .timeseries, .events, .stick_functions, .fits, .settings

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x));
p.addRequired('events_tsv', @(x) ischar(x) || isstring(x));
p.addParameter('TR', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('MaskNii', '', @(x) ischar(x) || isstring(x));
p.addParameter('Conditions', {}, @(x) ischar(x) || iscellstr(x) || isstring(x));
p.addParameter('WindowSeconds', 30, @(x) isscalar(x) && x > 0);
p.addParameter('Models', {'logit','fir','sfir','canonical','spline','nlgamma'}, @(x) iscell(x) || isstring(x));
p.addParameter('ModelDependencyPolicy', 'skip', @(x) ischar(x) || isstring(x));
p.addParameter('OutputMat', '', @(x) ischar(x) || isstring(x));
p.addParameter('OutputMatVersion', '-v7.3', @(x) ischar(x) || isstring(x));
p.addParameter('SaveWholeBrainInMat', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('SignalSource', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('ImageSet', 'all', @(x) ischar(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('SignatureName', '', @(x) ischar(x) || isstring(x));
p.addParameter('SignatureNames', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('MaxSignatures', inf, @(x) isscalar(x) && x >= 1);
p.addParameter('UseParallel', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('AtlasObj', [], @(x) isempty(x) || isa(x, 'atlas'));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MapNames', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('WriteWholeBrain', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('WholeBrainOutputPrefix', '', @(x) ischar(x) || isstring(x));
p.addParameter('WholeBrainMode', 'auto', @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('WholeBrainOverwrite', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('WholeBrainPThresh', [], @(x) isempty(x) || (isscalar(x) && x > 0 && x < 1));
p.addParameter('WholeBrainThreshType', 'unc', @(x) ischar(x) || isstring(x));
p.addParameter('WholeBrainWriteThresholdedT', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('WholeBrainChunkSize', 50000, @(x) isscalar(x) && x >= 1);
p.addParameter('WholeBrainScaleMode', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('ReuseWholeBrainOutput', false, @(x) islogical(x) || isnumeric(x));
p.parse(fmri_nii, events_tsv, varargin{:});
opts = p.Results;
write_wholebrain = logical(opts.WriteWholeBrain) || ~isempty(char(opts.WholeBrainOutputPrefix));
if write_wholebrain
    wholebrain_modes = local_resolve_wholebrain_modes(opts.WholeBrainMode, opts.Models);
else
    wholebrain_modes = {};
end

% Convenience: if a signature name is provided, force signature mode
if strcmpi(char(opts.SignalSource), 'mean') && ~isempty(char(opts.SignatureName))
    warning('SignatureName provided while SignalSource=''mean''. Switching SignalSource to ''signature''.');
    opts.SignalSource = 'signature';
end

fmri_nii = char(fmri_nii);
events_tsv = char(events_tsv);
if ~exist(fmri_nii, 'file'), error('fMRI file not found: %s', fmri_nii); end
if ~exist(events_tsv, 'file'), error('Events file not found: %s', events_tsv); end

[tr_from_hdr, n_tp] = local_get_tr_and_ntp(fmri_nii);
TR = opts.TR;
if isempty(TR), TR = tr_from_hdr; end
if isempty(TR) || TR <= 0
    error('Could not infer TR from NIfTI header. Pass ''TR'' explicitly.');
end

signal_source = lower(strtrim(char(opts.SignalSource)));
if strcmp(signal_source, 'signatures'), signal_source = 'signature'; end
if strcmp(signal_source, 'maps'), signal_source = 'imageset'; end
if strcmp(signal_source, 'roi'), signal_source = 'atlas'; end
signature_meta = struct('signal_source', signal_source);
fits_by_signature = struct();

if strcmpi(signal_source, 'mean')
    [tc, ~, ~] = hrf_extract_timeseries_from_nii(fmri_nii, char(opts.MaskNii));
elseif strcmpi(signal_source, 'signature')
    dat_for_maps = fmri_data(fmri_nii, 'noverbose');
    if isa(opts.ImageSet, 'image_vector')
        error('SignalSource=''signature'' expects ImageSet to be a named signature set. Use SignalSource=''imageset'' for an image_vector map set.');
    end
    if ~isempty(opts.SignatureName)
        % Fast/safe path for one signature to avoid loading unavailable images in image_set='all'
        image_set = char(opts.ImageSet);
        sig_name = char(opts.SignatureName);
        if strcmpi(image_set, 'all') && strcmpi(sig_name, 'NPS')
            image_set = 'npsplus';
        end
        [tc, one_meta] = hrf_extract_signature_timeseries(dat_for_maps, ...
            'SimilarityMetric', opts.SimilarityMetric, ...
            'ImageSet', image_set, ...
            'SignatureName', sig_name);
        all_tc = tc;
        signature_meta.available_signatures = {one_meta.selected_signature};
        signature_meta.selected_signature = one_meta.selected_signature;
        signature_meta.selected_signatures = {one_meta.selected_signature};
        signature_meta.n_signatures = 1;
        signature_meta.similarity_metric = one_meta.similarity_metric;
        signature_meta.image_set = one_meta.image_set;
    else
        [all_tc, signature_meta] = hrf_extract_all_signature_timeseries(dat_for_maps, ...
            'SimilarityMetric', opts.SimilarityMetric, ...
            'ImageSet', opts.ImageSet);

        sig_names = signature_meta.available_signatures;
        if ~isempty(opts.SignatureNames)
            req = cellstr(string(opts.SignatureNames));
            [tf, idx] = ismember(req, sig_names);
            sig_names = sig_names(idx(tf));
            all_tc = all_tc(:, idx(tf));
        end
        if isempty(sig_names)
            error('No signatures matched the requested SignatureNames.');
        end
        if isfinite(opts.MaxSignatures)
            n_keep = min(numel(sig_names), opts.MaxSignatures);
            sig_names = sig_names(1:n_keep);
            all_tc = all_tc(:, 1:n_keep);
        end
        tc = all_tc(:, 1);
        signature_meta.selected_signature = sig_names{1};
        signature_meta.selected_signatures = sig_names;
        signature_meta.n_signatures = numel(sig_names);
    end


elseif strcmpi(signal_source, 'atlas')
    if isempty(opts.AtlasObj)
        error('SignalSource=''atlas'' requires AtlasObj input (atlas object).');
    end
    [all_tc, signature_meta] = hrf_extract_roi_timeseries(fmri_nii, opts.AtlasObj, 'Regions', opts.Regions);
    sig_names = signature_meta.available_signatures;
    if isfinite(opts.MaxSignatures)
        n_keep = min(numel(sig_names), opts.MaxSignatures);
        sig_names = sig_names(1:n_keep);
        all_tc = all_tc(:, 1:n_keep);
    end
    tc = all_tc(:, 1);
    signature_meta.selected_signature = sig_names{1};
    signature_meta.selected_signatures = sig_names;
    signature_meta.n_signatures = numel(sig_names);
elseif strcmpi(signal_source, 'imageset')
    dat_for_maps = fmri_data(fmri_nii, 'noverbose');
    [all_tc, signature_meta] = hrf_extract_imageset_timeseries(dat_for_maps, opts.ImageSet, ...
        'MapNames', opts.MapNames, ...
        'SimilarityMetric', opts.SimilarityMetric);
    sig_names = signature_meta.available_signatures;
    if isempty(sig_names)
        error('No maps were available from the requested ImageSet/MapNames.');
    end
    if isfinite(opts.MaxSignatures)
        n_keep = min(numel(sig_names), opts.MaxSignatures);
        sig_names = sig_names(1:n_keep);
        all_tc = all_tc(:, 1:n_keep);
    end
    tc = all_tc(:, 1);
    signature_meta.selected_signature = sig_names{1};
    signature_meta.selected_signatures = sig_names;
    signature_meta.n_signatures = numel(sig_names);
else
    error('Unknown SignalSource: %s. Use ''mean'', ''signature'', ''imageset'', or ''atlas''.', char(opts.SignalSource));
end

E = hrf_load_events_tsv(events_tsv);
if isempty(opts.Conditions)
    cond_names = unique(E.trial_type, 'stable');
else
    cond_names = cellstr(string(opts.Conditions));
end

[Runc, condition_groups] = hrf_build_stick_functions(E, cond_names, TR, n_tp);
cond_names = {condition_groups.label};
if strcmpi(signal_source, 'signature') || strcmpi(signal_source, 'imageset') || strcmpi(signal_source, 'atlas')
    siglist = signature_meta.selected_signatures;
    nSig = numel(siglist);
    sigfields = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(siglist));
    signature_meta.selected_signature_fields = sigfields;
    fit_cells = cell(1, nSig);
    window_seconds = opts.WindowSeconds;
    models = opts.Models;
    model_dependency_policy = opts.ModelDependencyPolicy;
    use_parallel = logical(opts.UseParallel) && license('test', 'Distrib_Computing_Toolbox');

    if use_parallel
        if isempty(gcp('nocreate'))
            parpool;
        end
        parfor si = 1:nSig
            fit_cells{si} = hrf_fit_all_models(all_tc(:, si), TR, Runc, window_seconds, models, ...
                'DependencyPolicy', model_dependency_policy);
        end
    else
        for si = 1:nSig
            fit_cells{si} = hrf_fit_all_models(all_tc(:, si), TR, Runc, window_seconds, models, ...
                'DependencyPolicy', model_dependency_policy);
        end
    end

    for si = 1:nSig
        fits_by_signature.(sigfields{si}) = fit_cells{si};
    end
    selected_idx = find(strcmp(siglist, signature_meta.selected_signature), 1);
    fits = fit_cells{selected_idx};
else
    fits = hrf_fit_all_models(tc, TR, Runc, opts.WindowSeconds, opts.Models, ...
        'DependencyPolicy', opts.ModelDependencyPolicy);
    fits_by_signature.mean = fits;
    signature_meta.selected_signature = 'mean_bold';
    signature_meta.selected_signatures = {'mean_bold'};
    signature_meta.available_signatures = {'mean_bold'};
    signature_meta.n_signatures = 1;
end

results = struct();
results.timeseries = tc;
results.events = E;
results.conditions = cond_names;
results.condition_groups = condition_groups;
results.stick_functions = Runc;
results.fits = fits;
results.settings = struct('TR', TR, 'window_seconds', opts.WindowSeconds, ...
    'fmri_nii', fmri_nii, 'events_tsv', events_tsv, 'mask_nii', char(opts.MaskNii), ...
    'signal_source', char(opts.SignalSource), 'wholebrain_modes', {wholebrain_modes});
if isscalar(wholebrain_modes)
    results.settings.wholebrain_mode = wholebrain_modes{1};
end
results.signature_meta = signature_meta;
results.fits_by_signature = fits_by_signature;
if exist('all_tc', 'var')
    results.all_timeseries = all_tc;
    results.timeseries_by_signature = local_timeseries_by_signature(all_tc, signature_meta);
end

if write_wholebrain
    wholebrain_prefix = char(opts.WholeBrainOutputPrefix);
    if isempty(wholebrain_prefix)
        wholebrain_prefix = local_default_wholebrain_prefix(fmri_nii);
    end

    write_thresholded_t = logical(opts.WholeBrainWriteThresholdedT) || ~isempty(opts.WholeBrainPThresh);
    wholebrain_by_model = struct();
    for m = 1:numel(wholebrain_modes)
        mode_name = wholebrain_modes{m};
        model_prefix = local_wholebrain_model_prefix(wholebrain_prefix, mode_name, numel(wholebrain_modes));
        model_field = matlab.lang.makeValidName(lower(mode_name));
        if logical(opts.ReuseWholeBrainOutput) && local_has_wholebrain_outputs(model_prefix)
            wholebrain_by_model.(model_field) = local_load_wholebrain_outputs(model_prefix);
            fprintf('Reusing existing whole-brain HRF outputs: %s_*.nii\n', model_prefix);
        else
            wholebrain_by_model.(model_field) = hrf_fit_wholebrain_stats(fmri_nii, events_tsv, ...
                'TR', TR, ...
                'MaskNii', char(opts.MaskNii), ...
                'Conditions', cond_names, ...
                'WindowSeconds', opts.WindowSeconds, ...
                'Mode', mode_name, ...
                'OutputPrefix', model_prefix, ...
                'PThresh', opts.WholeBrainPThresh, ...
                'ThreshType', opts.WholeBrainThreshType, ...
                'WriteThresholdedT', write_thresholded_t, ...
                'ChunkSize', opts.WholeBrainChunkSize, ...
                'ScaleMode', opts.WholeBrainScaleMode, ...
                'Overwrite', logical(opts.WholeBrainOverwrite));
        end
    end
    results.wholebrain_by_model = wholebrain_by_model;
    if isscalar(wholebrain_modes)
        results.wholebrain = wholebrain_by_model.(matlab.lang.makeValidName(lower(wholebrain_modes{1})));
    end
end

if ~isempty(opts.OutputMat)
    local_save_results_mat(char(opts.OutputMat), results, opts);
end
end

function tf = local_has_wholebrain_outputs(prefix)
tf = exist([prefix '_beta.nii'], 'file') == 2 && exist([prefix '_t.nii'], 'file') == 2;
end

function wholebrain = local_load_wholebrain_outputs(prefix)
wholebrain = hrf_load_wholebrain_stats(prefix);
end

function modes = local_resolve_wholebrain_modes(mode_input, models)
if iscell(mode_input) || (isstring(mode_input) && numel(mode_input) > 1)
    requested = cellstr(string(mode_input));
else
    requested = {char(mode_input)};
end

if isscalar(requested) && strcmpi(requested{1}, 'auto')
    requested = cellstr(string(models));
end

supported = {'fir', 'sfir', 'canonical', 'spline'};
modes = {};
for i = 1:numel(requested)
    mode = lower(strtrim(char(requested{i})));
    if strcmp(mode, 'auto')
        nested = local_resolve_wholebrain_modes('auto', models);
        modes = [modes, nested]; %#ok<AGROW>
        continue
    end
    if ~ismember(mode, supported)
        if ismember(mode, {'logit', 'nlgamma'})
            warning('run_hrf_pipeline:SkippingWholeBrainModel', ...
                'Skipping whole-brain %s maps: nonlinear voxelwise model is not supported by the fast 4D writer.', mode);
            continue
        end
        error('Unknown WholeBrainMode/model: %s. Use auto, FIR, sFIR, canonical, or spline.', mode);
    end
    if strcmp(mode, 'sfir')
        modes{end + 1} = 'sFIR'; %#ok<AGROW>
    else
        modes{end + 1} = mode; %#ok<AGROW>
    end
end

modes = unique(modes, 'stable');
if isempty(modes)
    modes = {'FIR'};
end
end

function model_prefix = local_wholebrain_model_prefix(prefix, mode_name, n_modes)
if n_modes == 1
    model_prefix = prefix;
else
    model_prefix = sprintf('%s_%s', prefix, lower(char(mode_name)));
end
end

function local_save_results_mat(output_mat, results, opts)
save_version = char(opts.OutputMatVersion);
results_to_save = results;
if ~logical(opts.SaveWholeBrainInMat)
    if isfield(results_to_save, 'wholebrain')
        results_to_save.wholebrain_paths = results_to_save.wholebrain.paths;
        if isfield(results_to_save.wholebrain, 'metadata_table')
            results_to_save.wholebrain_metadata_table = results_to_save.wholebrain.metadata_table;
        end
        results_to_save = rmfield(results_to_save, 'wholebrain');
    end
    if isfield(results_to_save, 'wholebrain_by_model')
        fields = fieldnames(results_to_save.wholebrain_by_model);
        wholebrain_paths_by_model = struct();
        wholebrain_metadata_by_model = struct();
        for i = 1:numel(fields)
            wb = results_to_save.wholebrain_by_model.(fields{i});
            if isfield(wb, 'paths')
                wholebrain_paths_by_model.(fields{i}) = wb.paths;
            end
            if isfield(wb, 'metadata_table')
                wholebrain_metadata_by_model.(fields{i}) = wb.metadata_table;
            end
        end
        results_to_save.wholebrain_paths_by_model = wholebrain_paths_by_model;
        results_to_save.wholebrain_metadata_by_model = wholebrain_metadata_by_model;
        results_to_save = rmfield(results_to_save, 'wholebrain_by_model');
    end
end

results = results_to_save;
if isempty(save_version)
    save(output_mat, 'results');
else
    save(output_mat, 'results', save_version);
end
end

function [TR, n_tp] = local_get_tr_and_ntp(fmri_nii)
info = niftiinfo(fmri_nii);
if numel(info.ImageSize) < 4
    error('Expected 4D fMRI image, got %d dimensions.', numel(info.ImageSize));
end
n_tp = info.ImageSize(4);
TR = [];
if isfield(info, 'PixelDimensions') && numel(info.PixelDimensions) >= 4
    TR = info.PixelDimensions(4);
end
end

function prefix = local_default_wholebrain_prefix(fmri_nii)
[p, f, e] = fileparts(fmri_nii);
if strcmpi(e, '.gz')
    [~, f] = fileparts(f);
end
prefix = fullfile(p, [f '_hrf_wholebrain']);
end

function ts_by_signature = local_timeseries_by_signature(all_tc, signature_meta)
ts_by_signature = struct();
if ~isfield(signature_meta, 'selected_signatures') || isempty(signature_meta.selected_signatures)
    return
end

sig_names = cellstr(string(signature_meta.selected_signatures));
if isfield(signature_meta, 'selected_signature_fields') && ...
        numel(signature_meta.selected_signature_fields) >= numel(sig_names)
    sig_fields = cellstr(string(signature_meta.selected_signature_fields));
else
    sig_fields = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(sig_names));
end

n = min(numel(sig_fields), size(all_tc, 2));
for i = 1:n
    ts_by_signature.(sig_fields{i}) = all_tc(:, i);
end
end
