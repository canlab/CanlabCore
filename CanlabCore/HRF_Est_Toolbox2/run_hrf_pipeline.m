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
%   'Models'          : subset of {'logit','sfir','canonical','spline','nlgamma'}
%                       (default all)
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
p.addParameter('Conditions', {}, @(x) iscellstr(x) || isstring(x));
p.addParameter('WindowSeconds', 30, @(x) isscalar(x) && x > 0);
p.addParameter('Models', {'logit','sfir','canonical','spline','nlgamma'}, @iscell);
p.addParameter('OutputMat', '', @(x) ischar(x) || isstring(x));
p.addParameter('SignalSource', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('ImageSet', 'all', @(x) ischar(x) || isstring(x));
p.addParameter('SignatureName', '', @(x) ischar(x) || isstring(x));
p.addParameter('SignatureNames', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MaxSignatures', inf, @(x) isscalar(x) && x >= 1);
p.addParameter('UseParallel', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('MapNames', {}, @(x) iscell(x) || isstring(x));
p.addParameter('AtlasObj', [], @(x) isempty(x) || isa(x, 'atlas'));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x));
p.parse(fmri_nii, events_tsv, varargin{:});
opts = p.Results;

% Convenience: if a signature name is provided, force signature mode
if strcmpi(char(opts.SignalSource), 'mean') && ~isempty(char(opts.SignatureName))
    warning('SignatureName provided while SignalSource=''mean''. Switching SignalSource to ''signature''.');
    opts.SignalSource = 'signature';
end

fmri_nii = char(fmri_nii);
events_tsv = char(events_tsv);
if ~exist(fmri_nii, 'file'), error('fMRI file not found: %s', fmri_nii); end
if ~exist(events_tsv, 'file'), error('Events file not found: %s', events_tsv); end

[~, tr_from_hdr, n_tp] = hrf_extract_timeseries_from_nii(fmri_nii, char(opts.MaskNii));
TR = opts.TR;
if isempty(TR), TR = tr_from_hdr; end

signal_source = lower(strtrim(char(opts.SignalSource)));
if strcmp(signal_source, 'signatures'), signal_source = 'signature'; end
if strcmp(signal_source, 'maps'), signal_source = 'imageset'; end
if strcmp(signal_source, 'roi'), signal_source = 'atlas'; end
signature_meta = struct('signal_source', signal_source);
fits_by_signature = struct();

if strcmpi(signal_source, 'mean')
    [tc, ~, ~] = hrf_extract_timeseries_from_nii(fmri_nii, char(opts.MaskNii));
elseif strcmpi(signal_source, 'signature')
    if ~isempty(opts.SignatureName)
        % Fast/safe path for one signature to avoid loading unavailable images in image_set='all'
        image_set = char(opts.ImageSet);
        sig_name = char(opts.SignatureName);
        if strcmpi(image_set, 'all') && strcmpi(sig_name, 'NPS')
            image_set = 'npsplus';
        end
        [tc, one_meta] = hrf_extract_signature_timeseries(fmri_nii, ...
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
        [all_tc, signature_meta] = hrf_extract_all_signature_timeseries(fmri_nii, ...
            'SimilarityMetric', opts.SimilarityMetric, ...
            'ImageSet', opts.ImageSet);

        sig_names = signature_meta.available_signatures;
        if ~isempty(opts.SignatureNames)
            req = cellstr(string(opts.SignatureNames));
            [tf, idx] = ismember(req, sig_names);
            sig_names = sig_names(idx(tf));
            all_tc = all_tc(:, idx(tf));
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
    [all_tc, signature_meta] = hrf_extract_imageset_timeseries(fmri_nii, char(opts.ImageSet), ...
        'MapNames', opts.MapNames, ...
        'SimilarityMetric', opts.SimilarityMetric);
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
else
    error('Unknown SignalSource: %s. Use ''mean'', ''signature'', ''imageset'', or ''atlas''.', char(opts.SignalSource));
end

E = hrf_load_events_tsv(events_tsv);
if isempty(opts.Conditions)
    cond_names = unique(E.trial_type, 'stable');
else
    cond_names = cellstr(opts.Conditions);
end

Runc = hrf_build_stick_functions(E, cond_names, TR, n_tp);
if strcmpi(signal_source, 'signature') || strcmpi(signal_source, 'imageset') || strcmpi(signal_source, 'atlas')
    siglist = signature_meta.selected_signatures;
    nSig = numel(siglist);
    fit_cells = cell(1, nSig);
    use_parallel = logical(opts.UseParallel) && license('test', 'Distrib_Computing_Toolbox');

    if use_parallel
        if isempty(gcp('nocreate'))
            parpool;
        end
        parfor si = 1:nSig
            fit_cells{si} = hrf_fit_all_models(all_tc(:, si), TR, Runc, opts.WindowSeconds, opts.Models);
        end
    else
        for si = 1:nSig
            fit_cells{si} = hrf_fit_all_models(all_tc(:, si), TR, Runc, opts.WindowSeconds, opts.Models);
        end
    end

    for si = 1:nSig
        fits_by_signature.(siglist{si}) = fit_cells{si};
    end
    fits = fits_by_signature.(signature_meta.selected_signature);
else
    fits = hrf_fit_all_models(tc, TR, Runc, opts.WindowSeconds, opts.Models);
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
results.stick_functions = Runc;
results.fits = fits;
results.settings = struct('TR', TR, 'window_seconds', opts.WindowSeconds, ...
    'fmri_nii', fmri_nii, 'events_tsv', events_tsv, 'mask_nii', char(opts.MaskNii), ...
    'signal_source', char(opts.SignalSource));
results.signature_meta = signature_meta;
results.fits_by_signature = fits_by_signature;

if ~isempty(opts.OutputMat)
    save(char(opts.OutputMat), 'results');
end
end
