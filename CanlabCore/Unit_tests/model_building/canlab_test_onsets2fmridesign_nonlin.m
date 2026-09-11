function tests = canlab_test_onsets2fmridesign_nonlin
%CANLAB_TEST_ONSETS2FMRIDESIGN_NONLIN Numerical invariants for nonlinear VNL convolution.
%
% Covers onsets2fmridesign_nonlin and the maintenance fixes to
% modifiedconv_wager2005. The key anchors are:
%
%   - an isolated event has peak height 1
%   - a prior stimulus exactly 1 s earlier reproduces the original x = 2
%     saturation, so the time-discounted model agrees with the published
%     equations at the rate they were fit at
%   - with the nonlinearity switched off, the function reproduces
%     onsets2fmridesign bit for bit (validates neural-function construction
%     and downsampling)
%   - with 'condf', 'notimediscount' and res = 1/TR, it reproduces
%     modifiedconv_wager2005 bit for bit
%
% Requires SPM on the path (spm_hrf).

tests = functiontests(localfunctions);
end


function setupOnce(tc)
tc.assumeTrue(logical(exist('spm_hrf', 'file')), 'SPM (spm_hrf) is not on the path.');
end


% -------------------------------------------------------------------------
% Scaling and history index
% -------------------------------------------------------------------------

function test_isolated_event_peaks_at_one(tc)

X = onsets2fmridesign_nonlin({0}, 0.1, 40);
tc.verifyEqual(max(X(:, 1)), 1, 'AbsTol', 1e-3);

end


function test_one_second_lag_reproduces_x_equals_2(tc)
% The equations were fit to 1 Hz trains, so a single prior stimulus 1 s back
% must land exactly on x = 2.

[~, ~, ~, ~, out] = onsets2fmridesign_nonlin({[0; 1]}, 0.1, 40);
tc.verifyEqual(out.x_events{1}(2), 2, 'AbsTol', 1e-10);

% and an isolated event on x = 1
[~, ~, ~, ~, out] = onsets2fmridesign_nonlin({0}, 0.1, 40);
tc.verifyEqual(out.x_events{1}(1), 1, 'AbsTol', 1e-10);

end


function test_discount_function(tc)

[~, ~, ~, ~, out] = onsets2fmridesign_nonlin({0}, 1, 20);
w = out.discount_fcn;

tc.verifyEqual(w(1), 1, 'AbsTol', 1e-12);           % anchor at 1 s
tc.verifyEqual(w(0.25), 1, 'AbsTol', 1e-12);        % held at 1 below 1 s
tc.verifyEqual(w(5), 0.10, 'AbsTol', 1e-10);        % Miezin et al. 2000
tc.verifyLessThan(w(7), 0.05);                      % ~nothing left by the cutoff

% Beyond 'maxhistory', a prior stimulus has no effect at all
[~, ~, ~, ~, o8] = onsets2fmridesign_nonlin({[0; 8]}, 0.1, 40);
tc.verifyEqual(o8.x_events{1}(2), 1, 'AbsTol', 1e-12);

end


function test_history_index_is_monotonic_in_lag(tc)

lags = [1 2 3 4 5 6];
x = zeros(size(lags));

for i = 1:length(lags)
    [~, ~, ~, ~, out] = onsets2fmridesign_nonlin({[0; lags(i)]}, 0.1, 60);
    x(i) = out.x_events{1}(2);
end

tc.verifyTrue(all(diff(x) < 0), 'History index must decrease with longer lags.');
tc.verifyGreaterThanOrEqual(min(x), 1);

end


function test_hrf_fcn_matches_returned_hrf(tc)
% out.hrf_fcn(1) is the isolated-event HRF returned as the 4th output

[~, ~, ~, hrf, out] = onsets2fmridesign_nonlin({0}, 0.1, 40);

tc.verifyEqual(out.hrf_fcn(1), hrf, 'AbsTol', 1e-12);
tc.verifyEqual(max(hrf), 1, 'AbsTol', 1e-12);

% Higher history index means a lower, later-peaking response
h2 = out.hrf_fcn(2);
tc.verifyLessThan(max(h2), max(hrf));

[~, i1] = max(hrf);
[~, i2] = max(h2);
tc.verifyGreaterThan(i2, i1);

end


% -------------------------------------------------------------------------
% Agreement with the functions this generalizes
% -------------------------------------------------------------------------

function test_matches_onsets2fmridesign_when_nonlinearity_is_off(tc)
% Constant height, no shape modulation => plain linear convolution.

TR = 1; runlen = 120;
noneq = {'heighteq', @(x) ones(size(x)), 'delayeq', [], 'ttopeakeq', []};

% Impulse events
ons = {[5 20 40 70 95]'};
Xlin = onsets2fmridesign(ons, TR, runlen);
Xnl = onsets2fmridesign_nonlin(ons, TR, runlen, [], noneq{:});
tc.verifyEqual(Xnl, Xlin, 'AbsTol', 1e-12);

% Epochs (durations in column 2)
onsd = {[5 3; 20 3; 40 3; 70 3; 95 3]};
Xlin = onsets2fmridesign(onsd, TR, runlen);
Xnl = onsets2fmridesign_nonlin(onsd, TR, runlen, [], noneq{:});
tc.verifyEqual(Xnl, Xlin, 'AbsTol', 1e-12);

% Epochs via the 'durs' keyword
Xlin = onsets2fmridesign(ons, TR, runlen, [], 'durs', 3);
Xnl = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'durs', 3, noneq{:});
tc.verifyEqual(Xnl, Xlin, 'AbsTol', 1e-12);

end


function test_matches_modifiedconv_wager2005(tc)
% 'condf' input + count-based history + HRF sampled at the TR is exactly the
% original implementation.

TR = 1;
condf = [1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]';

Xold = modifiedconv_wager2005(TR, condf);
Xnew = onsets2fmridesign_nonlin(condf, TR, [], [], 'condf', 'notimediscount', 'res', 1 / TR);

tc.verifyEqual(Xnew(:, 1:end-1), Xold, 'AbsTol', 1e-12);
tc.verifyEqual(Xnew(:, end), ones(size(Xnew, 1), 1));   % intercept is appended

end


function test_all_three_equations_are_applied(tc)
% Height, onset delay, and time to peak must each change the prediction.

TR = 1;
condf = [ones(5, 1); zeros(15, 1)];

Xall = modifiedconv_wager2005(TR, condf);
Xnod = modifiedconv_wager2005(TR, condf, 'delay', []);
Xnop = modifiedconv_wager2005(TR, condf, 'ttopeak', []);
Xnoh = modifiedconv_wager2005(TR, condf, 'height', []);

tc.verifyGreaterThan(max(abs(Xall - Xnod)), 0.01, 'Onset delay equation has no effect.');
tc.verifyGreaterThan(max(abs(Xall - Xnop)), 0.01, 'Time-to-peak equation has no effect.');
tc.verifyGreaterThan(max(abs(Xall - Xnoh)), 0.01, 'Height equation has no effect.');

end


function test_modifiedconv_wager2005_fractional_tr(tc)
% Used to warn ("Integer operands are required for colon operator") and
% silently truncate the history window.

condf = [1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]';

X = tc.verifyWarningFree(@() modifiedconv_wager2005(1.3, condf));
tc.verifySize(X, [length(condf) 1]);

end


% -------------------------------------------------------------------------
% Multiple event types
% -------------------------------------------------------------------------

function test_no_cross_event_type_influence_by_default(tc)

TR = 1;
ons = {[5 40]' [8 45]'};    % type 2 always follows type 1 by 3 s

[~, ~, ~, ~, out] = onsets2fmridesign_nonlin(ons, TR, 120);

tc.verifyEqual(cat(1, out.x_events{:}), ones(4, 1), 'AbsTol', 1e-12);

end


function test_crossweights_enable_cross_type_influence(tc)

TR = 1;
ons = {[5 40]' [8 45]'};

[~, ~, ~, ~, out] = onsets2fmridesign_nonlin(ons, TR, 120, [], 'crossweights', [1 0.5; 0.5 1]);

tc.verifyEqual(out.x_events{1}, ones(2, 1), 'AbsTol', 1e-12);   % nothing precedes type 1
tc.verifyGreaterThan(min(out.x_events{2}), 1);                  % type 1 saturates type 2

end


% -------------------------------------------------------------------------
% Output shapes and onsets2fmridesign parity
% -------------------------------------------------------------------------

function test_output_shapes(tc)

TR = 1; runlen = 120;
ons = {[5 20 40 70 95]' [10 30 55 80]'};

Xlin = onsets2fmridesign(ons, TR, runlen);
[X, delta, delta_hires, hrf, out] = onsets2fmridesign_nonlin(ons, TR, runlen);

tc.verifySize(X, size(Xlin));                       % same shape, intercept included
tc.verifyEqual(X(:, end), ones(size(X, 1), 1));     % intercept is last
tc.verifyNumElements(delta, 2);
tc.verifyNumElements(delta{1}, size(X, 1));
tc.verifyEqual(sum(delta{1}), 5);
tc.verifySize(delta_hires, [runlen * 16, 2]);
tc.verifyEqual(max(hrf), 1, 'AbsTol', 1e-6);        % unit-history HRF
tc.verifyTrue(out.dotimediscount);

end


function test_saturates_relative_to_linear(tc)

TR = 1; runlen = 200;
ons = {(0:2:60)'};   % dense 2 s ISI train

Xlin = onsets2fmridesign(ons, TR, runlen);
Xnl = onsets2fmridesign_nonlin(ons, TR, runlen);

tc.verifyLessThan(max(Xnl(:, 1)), max(Xlin(:, 1)));

% Isolated events at the same amplitude scale in both
onsi = {[10 60 110]'};
Xlin = onsets2fmridesign(onsi, TR, runlen);
Xnl = onsets2fmridesign_nonlin(onsi, TR, runlen);
tc.verifyEqual(max(Xnl(:, 1)), max(Xlin(:, 1)), 'AbsTol', 1e-3);

end


function test_length_and_onset_handling_matches_onsets2fmridesign(tc)

TR = 2;
ons = {[10 40 70]'};

% Explicit length
X = onsets2fmridesign_nonlin(ons, TR, 120);
tc.verifySize(X, [60, 2]);

% No length: grows to cover the onsets, same rule as onsets2fmridesign
X = onsets2fmridesign_nonlin(ons, TR);
Xlin = onsets2fmridesign(ons, TR);
tc.verifySize(X, size(Xlin));

% Onsets past the end are dropped
X = onsets2fmridesign_nonlin({[10 40 700]'}, TR, 120);
Xref = onsets2fmridesign_nonlin({[10 40]'}, TR, 120);
tc.verifyEqual(X, Xref, 'AbsTol', 1e-12);

% Fractional TR
X = onsets2fmridesign_nonlin({[10 30]'}, 1.3, 130);
tc.verifySize(X, [100, 2]);

end


function test_durs_keyword_matches_column2(tc)

TR = 2; runlen = 120;
ons = {[10 40 70]' [20 55 95]'};
onsd = cellfun(@(x) [x repmat(4, size(x, 1), 1)], ons, 'UniformOutput', false);

tc.verifyEqual(onsets2fmridesign_nonlin(ons, TR, runlen, [], 'durs', 4), ...
    onsets2fmridesign_nonlin(onsd, TR, runlen), 'AbsTol', 1e-12);

% Per-trial durations via a cell array
d = {[2 4 6]' [1 3 5]'};
tc.verifyEqual(onsets2fmridesign_nonlin(ons, TR, runlen, [], 'durs', d), ...
    onsets2fmridesign_nonlin({[ons{1} d{1}] [ons{2} d{2}]}, TR, runlen), 'AbsTol', 1e-12);

% Numeric and row-vector onsets
tc.verifyEqual(onsets2fmridesign_nonlin([10 40 70], TR, runlen, [], 'durs', 4), ...
    onsets2fmridesign_nonlin({[10 40 70]'}, TR, runlen, [], 'durs', 4), 'AbsTol', 1e-12);

end


function test_keywords_without_positional_args(tc)

TR = 2;
ons = {[10 40 70]'};

tc.verifyEqual(onsets2fmridesign_nonlin(ons, TR, 'durs', 4), ...
    onsets2fmridesign_nonlin(ons, TR, [], [], 'durs', 4), 'AbsTol', 1e-12);

end


function test_rejects_incompatible_options(tc)

TR = 1; runlen = 120;
ons = {[5 20 40]'};

tc.verifyError(@() onsets2fmridesign_nonlin(ons, TR, runlen, spm_hrf(1)), '');
tc.verifyError(@() onsets2fmridesign_nonlin(ons, TR, runlen, [], 'parametric_standard', {[1 2 3]'}), '');
tc.verifyError(@() onsets2fmridesign_nonlin(ons, TR, runlen, [], 'customneural', ones(5, 1)), '');
tc.verifyError(@() onsets2fmridesign_nonlin(ons, TR, [10 20]), '');

end


function test_noundershoot(tc)

[~, ~, ~, h_default] = onsets2fmridesign_nonlin({0}, 1, 40);
[~, ~, ~, h_no_us] = onsets2fmridesign_nonlin({0}, 1, 40, [], 'noundershoot');

tc.verifyLessThan(min(h_default), -0.01);
tc.verifyGreaterThanOrEqual(min(h_no_us), -1e-12);

end


function test_nonorm(tc)
% Without normalization the isolated-event peak is heighteq(1)

heighteq = @(x) 1.7141 .* exp(-2.1038 .* x) + 0.4932 .* exp(-0.0770 .* x);

X = onsets2fmridesign_nonlin({10}, 0.1, 60, [], 'nonorm');
tc.verifyEqual(max(X(:, 1)), heighteq(1), 'AbsTol', 1e-3);

end


% -------------------------------------------------------------------------
% Single trial
% -------------------------------------------------------------------------

function test_singletrial_shape(tc)

TR = 1; runlen = 120;
ons = {[10 40 70 100]'};

X = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'singletrial');

tc.verifySize(X, [runlen / TR, 5]);   % 4 trials + intercept

% Well-separated trials each look like an isolated event
Xone = onsets2fmridesign_nonlin({10}, TR, runlen);
tc.verifyEqual(X(:, 1), Xone(:, 1), 'AbsTol', 1e-12);

end


function test_singletrial_history_accumulates_across_trials(tc)
% Splitting into one regressor per trial must NOT reset the stimulus
% history: a dense train's later trials are still saturated. (An earlier
% version computed history per regressor, so every trial got x = 1.)

TR = 1; runlen = 120;
ons = {(10:2:24)'};        % dense 2 s ISI train

[X, ~, ~, ~, out] = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'singletrial');

x_per_trial = cellfun(@(v) v(1), out.x_events);

tc.verifyEqual(x_per_trial(1), 1, 'AbsTol', 1e-12);         % nothing precedes trial 1
tc.verifyTrue(all(diff(x_per_trial) >= -1e-12), 'History must not decrease across trials.');
tc.verifyTrue(all(diff(x_per_trial(1:4)) > 0), 'History must build up across early trials.');
tc.verifyGreaterThan(x_per_trial(end), 1.5);

% Only stimuli within 'maxhistory' contribute, so x reaches a plateau rather
% than growing without limit. At a 2 s ISI that is 1 + w(2) + w(4) + w(6).
w = out.discount_fcn;
tc.verifyEqual(x_per_trial(end), 1 + w(2) + w(4) + w(6), 'AbsTol', 1e-10);

% Later trials therefore have smaller amplitudes than the first
tc.verifyLessThan(max(X(:, end-1)), max(X(:, 1)));

% And the single-trial regressors sum to the grouped model
Xgrouped = onsets2fmridesign_nonlin(ons, TR, runlen);
tc.verifyEqual(sum(X(:, 1:end-1), 2), Xgrouped(:, 1), 'AbsTol', 1e-12);

end


function test_singletrial_with_durations(tc)

TR = 1; runlen = 120;
ons = {[10 60]'};

X = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'singletrial', 'durs', 3);

tc.verifySize(X, [runlen / TR, 3]);   % 2 trials + intercept

Xone = onsets2fmridesign_nonlin({[10 3]}, TR, runlen);
tc.verifyEqual(X(:, 1), Xone(:, 1), 'AbsTol', 1e-12);

end


% -------------------------------------------------------------------------
% Custom neural drive via 'deltahires'
% -------------------------------------------------------------------------

function test_deltahires_input(tc)

TR = 1; runlen = 60; res = 16;

drive = zeros(runlen * res, 2);
drive(10 * res + 1, 1) = 1;
drive(30 * res + 1, 2) = 1;

X = onsets2fmridesign_nonlin(drive, TR, [], [], 'deltahires');
Xref = onsets2fmridesign_nonlin({10, 30}, TR, runlen);

tc.verifySize(X, [runlen / TR, 3]);
tc.verifyEqual(X, Xref, 'AbsTol', 1e-12);

% 'singletrial' has no trial structure to work with here
tc.verifyError(@() onsets2fmridesign_nonlin(drive, TR, [], [], 'deltahires', 'singletrial'), '');

end


% -------------------------------------------------------------------------
% Error paths
% -------------------------------------------------------------------------

function test_invalid_onsets_error(tc)

TR = 1; runlen = 60;

tc.verifyError(@() onsets2fmridesign_nonlin({[-5 10]'}, TR, runlen), '');
tc.verifyError(@() onsets2fmridesign_nonlin({[]}, TR, runlen), '');

end
