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
p.parse(fmri_nii, events_tsv, varargin{:});
opts = p.Results;

fmri_nii = char(fmri_nii);
events_tsv = char(events_tsv);
if ~exist(fmri_nii, 'file'), error('fMRI file not found: %s', fmri_nii); end
if ~exist(events_tsv, 'file'), error('Events file not found: %s', events_tsv); end

[tc, tr_from_hdr, n_tp] = hrf_extract_timeseries_from_nii(fmri_nii, char(opts.MaskNii));
TR = opts.TR;
if isempty(TR), TR = tr_from_hdr; end

E = hrf_load_events_tsv(events_tsv);
if isempty(opts.Conditions)
    cond_names = unique(E.trial_type, 'stable');
else
    cond_names = cellstr(opts.Conditions);
end

Runc = hrf_build_stick_functions(E, cond_names, TR, n_tp);
fits = hrf_fit_all_models(tc, TR, Runc, opts.WindowSeconds, opts.Models);

results = struct();
results.timeseries = tc;
results.events = E;
results.conditions = cond_names;
results.stick_functions = Runc;
results.fits = fits;
results.settings = struct('TR', TR, 'window_seconds', opts.WindowSeconds, ...
    'fmri_nii', fmri_nii, 'events_tsv', events_tsv, 'mask_nii', char(opts.MaskNii));

if ~isempty(opts.OutputMat)
    save(char(opts.OutputMat), 'results');
end
end
