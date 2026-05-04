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
p.addParameter('SimilarityMetric', 'dot_product', @(x) ischar(x) || isstring(x));
p.addParameter('ImageSet', 'all', @(x) ischar(x) || isstring(x));
p.addParameter('SignatureName', '', @(x) ischar(x) || isstring(x));
p.parse(fmri_nii, events_tsv, varargin{:});
opts = p.Results;

fmri_nii = char(fmri_nii);
events_tsv = char(events_tsv);
if ~exist(fmri_nii, 'file'), error('fMRI file not found: %s', fmri_nii); end
if ~exist(events_tsv, 'file'), error('Events file not found: %s', events_tsv); end

[~, tr_from_hdr, n_tp] = hrf_extract_timeseries_from_nii(fmri_nii, char(opts.MaskNii));
TR = opts.TR;
if isempty(TR), TR = tr_from_hdr; end

signal_source = lower(char(opts.SignalSource));
signature_meta = struct();
fits_by_signature = struct();

if strcmpi(signal_source, 'mean')
    [tc, ~, ~] = hrf_extract_timeseries_from_nii(fmri_nii, char(opts.MaskNii));
elseif strcmpi(signal_source, 'signature')
    [all_tc, signature_meta] = hrf_extract_all_signature_timeseries(fmri_nii, ...
        'SimilarityMetric', opts.SimilarityMetric, ...
        'ImageSet', opts.ImageSet);

    sig_names = signature_meta.available_signatures;
    if ~isempty(opts.SignatureName)
        sel = find(strcmp(sig_names, char(opts.SignatureName)), 1);
        if isempty(sel)
            error('Requested SignatureName "%s" not found.', char(opts.SignatureName));
        end
        sig_names = sig_names(sel);
        all_tc = all_tc(:, sel);
    end

    tc = all_tc(:, 1);
    signature_meta.selected_signature = sig_names{1};
    signature_meta.selected_signatures = sig_names;
    signature_meta.n_signatures = numel(sig_names);
else
    error('Unknown SignalSource: %s. Use ''mean'' or ''signature''.', char(opts.SignalSource));
end

E = hrf_load_events_tsv(events_tsv);
if isempty(opts.Conditions)
    cond_names = unique(E.trial_type, 'stable');
else
    cond_names = cellstr(opts.Conditions);
end

Runc = hrf_build_stick_functions(E, cond_names, TR, n_tp);
if strcmpi(signal_source, 'signature')
    for si = 1:numel(signature_meta.selected_signatures)
        sig = signature_meta.selected_signatures{si};
        sig_idx = find(strcmp(signature_meta.available_signatures, sig), 1);
        sig_tc = all_tc(:, sig_idx);
        fits_by_signature.(sig) = hrf_fit_all_models(sig_tc, TR, Runc, opts.WindowSeconds, opts.Models);
    end
    fits = fits_by_signature.(signature_meta.selected_signature);
else
    fits = hrf_fit_all_models(tc, TR, Runc, opts.WindowSeconds, opts.Models);
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
