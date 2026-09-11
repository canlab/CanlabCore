function tests = canlab_test_onsets2fmridesign
%CANLAB_TEST_ONSETS2FMRIDESIGN Option-by-option checks for onsets2fmridesign.
%
% Covers the major use cases: run length and TR handling, onset placement,
% epochs (both routes), the HRF/basis-set options, parametric modulators,
% single-trial models, custom neural functions, and the error paths.
%
% The 'durs' tests are regression tests for a 2026 fix. Before it, the
% keyword was unusable in every form: the cell branch replaced the caller's
% durations with ons ./ TR; scalar and per-condition durations were divided
% by the TR even though the epoch builder reads them as seconds; a plain
% numeric onset vector threw "Brace indexing is not supported"; and an
% unrecognized format returned empty onsets. The invariant those tests pin
% down is that 'durs' and column 2 of ons are two spellings of the same
% thing, both in seconds.
%
% Requires SPM on the path (spm_hrf, spm_get_bf).

tests = functiontests(localfunctions);
end


function setupOnce(tc)
tc.assumeTrue(logical(exist('spm_hrf', 'file')), 'SPM (spm_hrf) is not on the path.');
end


% -------------------------------------------------------------------------
% Shapes, timing, and run length
% -------------------------------------------------------------------------

function test_basic_shape_and_intercept(tc)

TR = 1; runlen = 60;
ons = {[10 30]' [20 45]'};

X = onsets2fmridesign(ons, TR, runlen);

tc.verifySize(X, [runlen / TR, 3]);                        % 2 conditions + intercept
tc.verifyEqual(X(:, end), ones(runlen / TR, 1), 'Last column must be the intercept.');

end


function test_onset_placement(tc)
% The event at 10 s must land in the 11th TR bin (time 0 is element 1) and
% at high-resolution sample 1 + 10*16.

TR = 1; runlen = 60; res = 16;
[~, delta, delta_hires] = onsets2fmridesign({[10 30]'}, TR, runlen);

tc.verifyNumElements(delta, 1);
tc.verifyNumElements(delta{1}, runlen / TR);
tc.verifyEqual(find(delta{1}), [10 / TR + 1; 30 / TR + 1]);

tc.verifySize(delta_hires, [runlen * res, 1]);
tc.verifyEqual(find(delta_hires), [10 * res + 1; 30 * res + 1]);
tc.verifyEqual(sum(delta_hires), 2);                       % unit drive per impulse

end


function test_impulse_peak_amplitude_and_timing(tc)
% An isolated impulse peaks at 1 (default amplitude scaling), roughly 5 s
% after onset for the canonical HRF.

TR = 1; runlen = 60;
X = onsets2fmridesign({10}, TR, runlen);

[pk, pki] = max(X(:, 1));
tc.verifyEqual(pk, 1, 'AbsTol', 1e-3);
tc.verifyEqual((pki - 1) * TR, 15, 'AbsTol', 1);           % onset 10 s + ~5 s to peak

end


function test_length_handling(tc)

TR = 2; ons = {[10 40 70]'};

% Explicit length
X = onsets2fmridesign(ons, TR, 120);
tc.verifySize(X, [60, 2]);

% No length: grows to cover the onsets, rounded up to a whole number of TRs
X = onsets2fmridesign(ons, TR);
tc.verifySize(X, [70 / TR, 2]);

end


function test_fractional_tr(tc)

TR = 1.3; runlen = 130;
[X, delta] = onsets2fmridesign({[10 30]'}, TR, runlen);

tc.verifySize(X, [runlen / TR, 2]);
tc.verifyNumElements(delta{1}, runlen / TR);

% A length that is not a whole number of TRs is an error, not a silent trim
tc.verifyError(@() onsets2fmridesign({10}, 1.3, 100), '');

end


function test_onsets_past_end_are_dropped(tc)

TR = 1; runlen = 60;
X = onsets2fmridesign({[10 30 200]'}, TR, runlen);
Xref = onsets2fmridesign({[10 30]'}, TR, runlen);

tc.verifySize(X, [runlen / TR, 2]);
tc.verifyEqual(X, Xref, 'AbsTol', 1e-12);

end


function test_conditions_are_built_independently(tc)

TR = 1; runlen = 60;
ons = {[10 30]' [20 45]'};

X = onsets2fmridesign(ons, TR, runlen);
X1 = onsets2fmridesign(ons(1), TR, runlen);
X2 = onsets2fmridesign(ons(2), TR, runlen);

tc.verifyEqual(X(:, 1), X1(:, 1), 'AbsTol', 1e-12);
tc.verifyEqual(X(:, 2), X2(:, 1), 'AbsTol', 1e-12);

end


% -------------------------------------------------------------------------
% Durations: 'durs' keyword vs column 2 of ons
% -------------------------------------------------------------------------

function test_durs_scalar_matches_column2(tc)

TR = 2; runlen = 120;
ons = {[10 40 70]' [20 55 95]'};
onsd = cellfun(@(x) [x repmat(4, size(x, 1), 1)], ons, 'UniformOutput', false);

X_keyword = onsets2fmridesign(ons, TR, runlen, [], 'durs', 4);
X_column = onsets2fmridesign(onsd, TR, runlen);

tc.verifyEqual(X_keyword, X_column, 'AbsTol', 1e-12);

end


function test_durs_cell_per_trial_matches_column2(tc)

TR = 2; runlen = 120;
ons = {[10 40 70]' [20 55 95]'};
d = {[2 4 6]' [1 3 5]'};
onsd = {[ons{1} d{1}] [ons{2} d{2}]};

X_keyword = onsets2fmridesign(ons, TR, runlen, [], 'durs', d);
X_column = onsets2fmridesign(onsd, TR, runlen);

tc.verifyEqual(X_keyword, X_column, 'AbsTol', 1e-12);

% Row-vector durations inside the cells must work too
X_row = onsets2fmridesign(ons, TR, runlen, [], 'durs', {[2 4 6] [1 3 5]});
tc.verifyEqual(X_row, X_column, 'AbsTol', 1e-12);

end


function test_durs_vector_per_condition_matches_column2(tc)

TR = 2; runlen = 120;
ons = {[10 40 70]' [20 55 95]'};
onsd = {[ons{1} repmat(2, 3, 1)] [ons{2} repmat(6, 3, 1)]};

X_keyword = onsets2fmridesign(ons, TR, runlen, [], 'durs', [2 6]);
X_column = onsets2fmridesign(onsd, TR, runlen);

tc.verifyEqual(X_keyword, X_column, 'AbsTol', 1e-12);

end


function test_durs_are_in_seconds_not_trs(tc)
% Durations are in seconds, so the underlying neural function must not
% depend on the TR. This is the sharpest form of the regression: the old
% code divided 'durs' by the TR, so these two differed by a factor of 2.

runlen = 120; res = 16; dur = 4;

[~, ~, dh_tr1] = onsets2fmridesign({[10 40]'}, 1, runlen, [], 'durs', dur);
[~, ~, dh_tr2] = onsets2fmridesign({[10 40]'}, 2, runlen, [], 'durs', dur);

tc.verifyEqual(dh_tr1, dh_tr2, 'AbsTol', 1e-12);

% And the epoch really is 'dur' seconds of neural drive at 1/res per sample
epoch = dh_tr1(10 * res + 1 : 10 * res + 1 + dur * res);
tc.verifyEqual(epoch, repmat(1 / res, dur * res + 1, 1), 'AbsTol', 1e-12);
tc.verifyEqual(nnz(dh_tr1), 2 * (dur * res + 1));

end


function test_durs_with_numeric_and_row_vector_onsets(tc)
% Plain numeric onsets used to throw "Brace indexing is not supported for
% variables of type double"; row vectors were appended in the wrong
% orientation.

TR = 2; runlen = 120;

X_cellcol = onsets2fmridesign({[10 40 70]'}, TR, runlen, [], 'durs', 4);
X_numcol = onsets2fmridesign([10 40 70]', TR, runlen, [], 'durs', 4);
X_numrow = onsets2fmridesign([10 40 70], TR, runlen, [], 'durs', 4);

tc.verifyEqual(X_numcol, X_cellcol, 'AbsTol', 1e-12);
tc.verifyEqual(X_numrow, X_cellcol, 'AbsTol', 1e-12);

end


function test_epoch_amplitude_scaling(tc)
% Longer epochs deliver more neural drive and produce larger responses;
% a ~1 s epoch is close to the unit single-event response.

TR = 1; runlen = 80;
durs = [1 2 4 8];
pk = zeros(size(durs));

for i = 1:length(durs)
    X = onsets2fmridesign({[10 durs(i)]}, TR, runlen);
    pk(i) = max(X(:, 1));
end

tc.verifyTrue(all(diff(pk) > 0), 'Longer epochs must produce larger responses.');
tc.verifyEqual(pk(1), 1, 'AbsTol', 0.1);            % 1 s epoch ~ unit response
tc.verifyLessThan(pk(4), 8);                        % strongly sublinear in duration

end


function test_durs_ignored_when_ons_already_has_durations(tc)

TR = 2; runlen = 120;
onsd = {[10 4; 40 4; 70 4]};

lastwarn('');
X = onsets2fmridesign(onsd, TR, runlen, [], 'durs', 99);
tc.verifySubstring(lastwarn, 'already contains durations');

tc.verifyEqual(X, onsets2fmridesign(onsd, TR, runlen), 'AbsTol', 1e-12);

end


function test_durs_size_mismatch_errors(tc)

TR = 1; runlen = 60;

% Wrong number of cells
tc.verifyError(@() onsets2fmridesign({[10 30]' [5 15]'}, TR, runlen, [], 'durs', {[1 2]'}), '');

% Wrong number of durations within a cell
tc.verifyError(@() onsets2fmridesign({[10 30]'}, TR, runlen, [], 'durs', {[1 2 3]'}), '');

end


function test_unrecognized_durs_does_not_destroy_onsets(tc)
% An unparseable 'durs' should warn and fall back to impulse events, not
% return an empty design.

TR = 1; runlen = 60;
ons = {[10 30 50]'};

lastwarn('');
X = onsets2fmridesign(ons, TR, runlen, [], 'durs', [1 2 3 4 5]);
tc.verifySubstring(lastwarn, 'unrecognized');

tc.verifyEqual(X, onsets2fmridesign(ons, TR, runlen), 'AbsTol', 1e-12);

end


% -------------------------------------------------------------------------
% Argument-position handling
% -------------------------------------------------------------------------

function test_keywords_without_positional_args(tc)
% onsets2fmridesign(ons, TR, 'durs', 4) used to assign 'durs' to len.

TR = 2;
ons = {[10 40 70]'};

X = onsets2fmridesign(ons, TR, 'durs', 4);
Xref = onsets2fmridesign(ons, TR, [], [], 'durs', 4);

tc.verifyEqual(X, Xref, 'AbsTol', 1e-12);

% and it really applied the durations, rather than falling back to impulses
tc.verifyGreaterThan(max(X(:, 1)), 1.5);

end


% -------------------------------------------------------------------------
% HRF and basis-set options
% -------------------------------------------------------------------------

function test_default_hrf_is_amplitude_scaled(tc)

TR = 1; runlen = 60; res = 16;
[~, ~, ~, hrf] = onsets2fmridesign({10}, TR, runlen);

tc.verifyNumElements(hrf, 32 * res + 1);        % 32 s kernel at 16 Hz
tc.verifyEqual(max(hrf), 1, 'AbsTol', 1e-12);
tc.verifyLessThan(min(hrf), 0);                 % canonical HRF has an undershoot

end


function test_noampscale(tc)

TR = 1; runlen = 60;
[X, ~, ~, hrf] = onsets2fmridesign({10}, TR, runlen, [], 'noampscale');

tc.verifyEqual(sum(hrf), 1, 'AbsTol', 1e-6);    % spm_hrf sums to 1, unscaled
tc.verifyLessThan(max(X(:, 1)), 0.1);

end


function test_noundershoot(tc)

TR = 1; runlen = 60;
[~, ~, ~, h_default] = onsets2fmridesign({10}, TR, runlen);
[~, ~, ~, h_no_us] = onsets2fmridesign({10}, TR, runlen, [], 'noundershoot');

tc.verifyLessThan(min(h_default), -0.01);
tc.verifyGreaterThanOrEqual(min(h_no_us), -1e-12);
tc.verifyEqual(max(h_no_us), 1, 'AbsTol', 1e-12);

end


function test_custom_hrf(tc)

TR = 1; runlen = 60; res = 16;

[X, ~, ~, hrf] = onsets2fmridesign({10}, TR, runlen, spm_hrf(1));

tc.verifyNumElements(hrf, 32 * res + 1);        % resampled from 1 s to 16 Hz
tc.verifyEqual(sum(hrf), 1, 'AbsTol', 1e-6);    % custom HRFs are scaled to sum 1
tc.verifySize(X, [runlen / TR, 2]);

% Peak of the response is still at the expected latency
[~, pki] = max(X(:, 1));
tc.verifyEqual((pki - 1) * TR, 15, 'AbsTol', 1.5);

end


function test_basis_sets(tc)

TR = 1; runlen = 60;
ons = {[10 30]' [20 45]'};

X = onsets2fmridesign(ons, TR, runlen, 'hrf (with time derivative)');
tc.verifySize(X, [runlen / TR, 2 * 2 + 1]);

X = onsets2fmridesign(ons, TR, runlen, 'hrf (with time and dispersion derivatives)');
tc.verifySize(X, [runlen / TR, 2 * 3 + 1]);

X = onsets2fmridesign(ons, TR, runlen, 'Finite Impulse Response');
tc.verifySize(X, [runlen / TR, 2 * 15 + 1]);    % order 15

end


function test_norm_orthogonalizes_basis_set(tc)

TR = 1; runlen = 60;
ons = {[10 30 50]'};

[~, ~, ~, hrf] = onsets2fmridesign(ons, TR, runlen, 'hrf (with time derivative)', 'norm');

tc.verifySize(hrf, [size(hrf, 1), 2]);
tc.verifyEqual(mean(hrf), [0 0], 'AbsTol', 1e-10);              % mean-centered
tc.verifyEqual(vecnorm(hrf), [1 1], 'AbsTol', 1e-10);           % unit L2 norm
tc.verifyEqual(hrf(:, 1)' * hrf(:, 2), 0, 'AbsTol', 1e-10);     % orthogonal

end


% -------------------------------------------------------------------------
% Model-structure options
% -------------------------------------------------------------------------

function test_singletrial(tc)

TR = 1; runlen = 60;
ons = {[10 30]' [20 45]'};

X = onsets2fmridesign(ons, TR, runlen, [], 'singletrial');

tc.verifySize(X, [runlen / TR, 5]);             % 4 trials + intercept

% Each single-trial regressor is the response to one event alone
Xone = onsets2fmridesign({10}, TR, runlen);
tc.verifyEqual(X(:, 1), Xone(:, 1), 'AbsTol', 1e-12);

end


function test_singletrial_with_durations(tc)

TR = 1; runlen = 60;
ons = {[10 30]'};

X = onsets2fmridesign(ons, TR, runlen, [], 'singletrial', 'durs', 3);

tc.verifySize(X, [runlen / TR, 3]);             % 2 trials + intercept

Xone = onsets2fmridesign({[10 3]}, TR, runlen);
tc.verifyEqual(X(:, 1), Xone(:, 1), 'AbsTol', 1e-12);

end


function test_parametric_standard(tc)
% Two regressors per event type: an unmodulated response, and one driven by
% the mean-centered modulator values.
%
% Note on scaling: getPredictors/pmconv builds the modulator regressor from
% wh:wh_end, which spans two high-resolution samples for an impulse event,
% so the modulator regressor carries about twice the amplitude of the
% unmodulated one. Regressor scaling does not affect model fit or t values,
% so the tests below check shape and the mean-centering invariants rather
% than absolute amplitude.

TR = 1; runlen = 60;
ons = {[10 30]' [20 45]'};
pm = {[1 2]' [2 1]'};

X = onsets2fmridesign(ons, TR, runlen, [], 'parametric_standard', pm);

tc.verifySize(X, [runlen / TR, 5]);             % 2 per condition + intercept

% First regressor of each pair is exactly the unmodulated response
Xplain = onsets2fmridesign(ons, TR, runlen);
tc.verifyEqual(X(:, 1), Xplain(:, 1), 'AbsTol', 1e-12);
tc.verifyEqual(X(:, 3), Xplain(:, 2), 'AbsTol', 1e-12);

% Modulator values are mean-centered: adding a constant to every value must
% leave both regressors untouched
Xshift = onsets2fmridesign(ons, TR, runlen, [], 'parametric_standard', ...
    cellfun(@(v) v + 10, pm, 'UniformOutput', false));
tc.verifyEqual(Xshift, X, 'AbsTol', 1e-12);

% The modulator regressor is linear in the modulator spread
Xdouble = onsets2fmridesign(ons, TR, runlen, [], 'parametric_standard', ...
    cellfun(@(v) 2 * v, pm, 'UniformOutput', false));
tc.verifyEqual(Xdouble(:, 2), 2 * X(:, 2), 'AbsTol', 1e-10);
tc.verifyEqual(Xdouble(:, 1), X(:, 1), 'AbsTol', 1e-12);

% With modulators [1 2], the regressor contrasts the second event against
% the first, so it tracks resp(30) - resp(10)
Xa = onsets2fmridesign({10}, TR, runlen);
Xb = onsets2fmridesign({30}, TR, runlen);
tc.verifyGreaterThan(corr(X(:, 2), Xb(:, 1) - Xa(:, 1)), 0.999);

end


function test_parametric_singleregressor(tc)

TR = 1; runlen = 60;
ons = {[10 30]'};

X1 = onsets2fmridesign(ons, TR, runlen, [], 'parametric_singleregressor', {[1 1]'});
X2 = onsets2fmridesign(ons, TR, runlen, [], 'parametric_singleregressor', {[2 2]'});

tc.verifySize(X1, [runlen / TR, 2]);            % one regressor per condition

% Doubling every modulator doubles the regressor
tc.verifyEqual(X2(:, 1), 2 * X1(:, 1), 'AbsTol', 1e-10);

end


function test_nonlinsaturation(tc)

TR = 1; runlen = 80;
ons = {(0:2:40)'};                              % dense train

Xlin = onsets2fmridesign(ons, TR, runlen);
Xsat = onsets2fmridesign(ons, TR, runlen, [], 'nonlinsaturation');

tc.verifyLessThan(max(Xsat(:, 1)), max(Xlin(:, 1)));

% 'nononlin' is the explicit off switch
Xoff = onsets2fmridesign(ons, TR, runlen, [], 'nononlin');
tc.verifyEqual(Xoff, Xlin, 'AbsTol', 1e-12);

end


function test_customneural(tc)

TR = 1; runlen = 60;
neural = [0 .5 1 .5 0]';

X = onsets2fmridesign({[10 30]'}, TR, runlen, [], 'customneural', neural);

tc.verifySize(X, [runlen / TR, 2]);

% Sustained neural drive produces a larger response than a bare impulse
Ximp = onsets2fmridesign({[10 30]'}, TR, runlen);
tc.verifyGreaterThan(max(X(:, 1)), max(Ximp(:, 1)));

% Custom neural functions cannot be combined with epochs
tc.verifyError(@() onsets2fmridesign({[10 3; 30 3]}, TR, runlen, [], 'customneural', neural), '');

end


% -------------------------------------------------------------------------
% plotDesign: 'durs' must mean the same thing there
% -------------------------------------------------------------------------

function test_plotdesign_durs_matches_onsets2fmridesign(tc)
% plotDesign carried its own copy of the 'durs' ./ TR bug, so its design and
% its drawn epoch boxes were TR times too short.

tc.assumeTrue(logical(exist('plotDesign', 'file')), 'plotDesign is not on the path.');

TR = 2;
ons = {[10 40 70]'};
onsd = {[10 4; 40 4; 70 4]};

fh = figure('Visible', 'off');
c = onCleanup(@() close(fh));

X_keyword = plotDesign(ons, [], TR, 'durs', 4, 'samefig');
X_column = plotDesign(onsd, [], TR, 'samefig');

tc.verifyEqual(X_keyword, X_column, 'AbsTol', 1e-12);

% ...and the same design onsets2fmridesign would build for 4 s epochs
len = size(X_keyword, 1) * TR;
X_o2f = onsets2fmridesign(onsd, TR, len);
tc.verifyEqual(X_keyword(:, 1), X_o2f(:, 1), 'AbsTol', 1e-12);

end


function test_plotdesign_accepts_cell_durs(tc)
% A cell array of durations used to error on the './ TR' line.

tc.assumeTrue(logical(exist('plotDesign', 'file')), 'plotDesign is not on the path.');

TR = 2;
ons = {[10 40 70]'};
d = {[2 4 6]'};

fh = figure('Visible', 'off');
c = onCleanup(@() close(fh));

X_keyword = plotDesign(ons, [], TR, 'durs', d, 'samefig');
X_column = plotDesign({[ons{1} d{1}]}, [], TR, 'samefig');

tc.verifyEqual(X_keyword, X_column, 'AbsTol', 1e-12);

end


% -------------------------------------------------------------------------
% Error and warning paths
% -------------------------------------------------------------------------

function test_invalid_onsets_error(tc)

TR = 1; runlen = 60;

tc.verifyError(@() onsets2fmridesign({[-5 10]'}, TR, runlen), '');   % negative onsets
tc.verifyError(@() onsets2fmridesign({[]}, TR, runlen), '');         % empty condition

end


function test_unknown_keyword_warns(tc)

TR = 1; runlen = 60;

lastwarn('');
onsets2fmridesign({10}, TR, runlen, [], 'bogus_option');
tc.verifySubstring(lastwarn, 'Unknown input string option');

end
