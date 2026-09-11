function [X, delta, delta_hires, hrf, out] = onsets2fmridesign_nonlin(ons, TR, varargin)
% Build a design matrix from onsets, using history-dependent nonlinear
% convolution (Wager et al. 2005 VNL equations with time discounting)
%
% :Usage:
% ::
%
%     [X, delta, delta_hires, hrf, out] = onsets2fmridesign_nonlin(onsets, TR, [len], [], [keywords])
%
% A drop-in variant of onsets2fmridesign.m that replaces linear convolution
% with a nonlinear model: each event is convolved with its own hemodynamic
% response function (HRF), whose height, onset delay, and time to peak
% depend on the stimulus history preceding it. Onset/duration formats, run
% length, single-trial expansion, output shapes, and the appended intercept
% column all match onsets2fmridesign.
%
% This generalizes modifiedconv_wager2005.m in three ways:
%
%   1. Onsets are given in seconds (with optional durations) rather than as
%      an indicator vector sampled at the TR. Neural events are built at
%      high resolution ('res', default 16 Hz) and downsampled to the TR.
%
%   2. Multiple event types are handled. History is computed *within* event
%      type: by default event type A does not saturate event type B. The
%      'crossweights' option exposes the machinery for cross-type
%      influences, but the default (identity) turns them off.
%
%   3. Nonlinearity is discounted by how long ago each preceding stimulus
%      occurred (on by default; see "Time-based discounting" below).
%
% :Nonlinearity model:
%
% The Wager et al. 2005 equations map a stimulus-history index x onto three
% HRF parameters. x = 1 is an isolated event; x = 2 is the second event of a
% 1 Hz train; and so on:
%
%     height(x)  =   1.7141*exp(-2.1038*x) + 0.4932*exp(-0.0770*x)
%     delay(x)   = -13.4097*exp(-1.0746*x) + 4.8733*exp(-0.1979*x)
%     ttopeak(x) =  37.5445*exp(-2.6760*x) - 3.2046*exp(-0.2120*x) + 5.6344
%
% height scales the peak amplitude, delay is spm_hrf parameter 6 (onset
% delay, s), and ttopeak is spm_hrf parameter 1 (time to peak, s).
% Predictors are rescaled by 1/height(1) so an isolated event peaks at 1.
%
% :Time-based discounting (default; 'notimediscount' to turn off):
%
% The equations above were estimated from 1 Hz trains, so x is only defined
% as an event count. To generalize to arbitrary onset times, each preceding
% stimulus contributes to x in proportion to an exponential decay in the
% time since it occurred:
%
%     x(t) = 1 + sum_k  s(t_k) * w(t - t_k),   w(dt) = exp(-lambda*(dt - 1))
%
% where s is the neural drive (1 per impulse event; 1/res per sample of an
% epoch, so a 1 s epoch contributes ~1 unit, matching onsets2fmridesign).
% w(1) = 1, so a single stimulus 1 s earlier reproduces the original x = 2
% saturation exactly. lambda is set from 'decay5' so that w(5) = 0.10, i.e.
% about 10% residual saturation from an event 5 s ago, consistent with the
% Miezin et al. (2000) review. w is held at 1 for dt < 1 s (the equations
% are not identified below the 1 Hz rate at which they were fit), and
% stimuli more than 'maxhistory' seconds back (default 7 s, where w = 0.03)
% are ignored.
%
% Because the discount kernel sums to a finite value, a sustained 1 Hz train
% drives x toward about 3.3 rather than counting up without limit, so peak
% height saturates near 58% of an isolated event rather than 34%. Use
% 'decay5' or 'lambda' to retune this.
%
% With 'notimediscount', x reverts to the original count-based index used by
% modifiedconv_wager2005: the number of same-type events in the preceding
% 12 s, scaled by their density in TR-sized bins. That path assumes impulse
% events sampled at the TR and is provided for comparison only.
%
% :Differences from onsets2fmridesign:
%
%   - No custom HRF or basis set. The nonlinear model is parameterized
%     through spm_hrf, so the 4th positional argument (HRF name or values)
%     must be empty. Basis-set options ('norm', 'noampscale') do not apply.
%   - No FIR / deconvolution matrix output. The FIR (DX) model assumes a
%     shape-invariant, linear response and is not compatible with a model in
%     which the response shape changes trial to trial; use
%     tor_make_deconv_mtx3 directly if you need one.
%   - No parametric modulators. Custom neural response functions are
%     supported through the 'deltahires' input mode rather than
%     'customneural'.
%   - 'nonlinsaturation' / 'nononlin' are ignored: nonlinearity is intrinsic
%     here. Use 'notimediscount' for the original count-based history index,
%     or onsets2fmridesign for a fully linear model.
%
% :Inputs:
%
%   **onsets:**
%        - 1st column is onsets for events in seconds
%        - 2nd column is optional durations for each event, in seconds
%        - Enter a single condition, or a cell vector with one cell per
%          event type.
%        With the 'condf' keyword, ons is instead an indicator vector
%        sampled at the TR (zeros and integer event codes), as accepted by
%        modifiedconv.m / modifiedconv_wager2005.m.
%        With the 'deltahires' keyword, ons is a [samples x n event types]
%        neural drive function already sampled at 'res' Hz.
%
%   **TR:**
%        Repetition time in seconds
%
% :Optional Inputs: First two are positional, then keywords:
%
%   **len:**
%        Length in sec for the model, usually n_images * TR, or [] to skip.
%
%   **(HRF):**
%        Must be [] -- present only so that calls written for
%        onsets2fmridesign work unchanged.
%
%   **'durs':**
%        Durations in seconds: a scalar (same for all events), one value per
%        event type, or a cell array with one duration per trial. Durations
%        may also be given as column 2 of ons.
%
%   **'noundershoot':**
%        Remove the post-stimulus undershoot (sets spm_hrf parameter 5, the
%        response-to-undershoot ratio, to Inf).
%
%   **'singletrial':**
%        One regressor per trial. History still accumulates across all
%        trials of the original event types, so each trial's shape reflects
%        the events that preceded it.
%
%   **'notimediscount':**
%        Use the original count-based stimulus-history index instead of
%        exponential time discounting. Aliases: 'original', 'nodiscount'.
%
%   **'decay5':**
%        Residual fraction of the 1 s nonlinear effect remaining for a
%        stimulus 5 s earlier. Default 0.10. Sets lambda = -log(decay5)/4.
%
%   **'lambda':**
%        Set the exponential decay rate directly, overriding 'decay5'.
%
%   **'maxhistory':**
%        Ignore stimuli more than this many seconds back. Default 7.
%
%   **'res':**
%        Samples per second for the internal neural/HRF representation.
%        Default 16, matching onsets2fmridesign.
%
%   **'crossweights':**
%        [n event types x n event types] matrix; element (i, j) weights the
%        history of event type j when computing x for event type i.
%        Default eye(k): no cross-type influence. Experimental.
%
%   **'heighteq', 'delayeq', 'ttopeakeq', 'uonseteq':**
%        Function handles overriding the equations above (uonseteq sets
%        spm_hrf parameter 2, undershoot delay; off by default). Enter [] to
%        disable modulation of that parameter.
%
%   **'hrfparams':**
%        Base spm_hrf parameter vector, [ttopeak uonset dispersion
%        udispersion rtou onsetdelay klength]. Default [6 16 1 1 6 0 32].
%
%   **'nonorm':**
%        Do not rescale by 1/height(1); predictors then carry the absolute
%        amplitudes predicted by the equations.
%
%   **'plot':**
%        Plot the resulting predictors.
%
%   **'verbose':**
%        Print a summary of the history index and decay parameters.
%
% :Outputs:
%
%   **X:**
%        Model sampled at TR, with an intercept appended as the last column.
%        Use X(:, 1:end-1) for the predictors alone.
%
%   **delta:**
%        Cell array of indicator vectors, sampled at TR (empty with
%        'deltahires', where there is no onset list to indicate)
%
%   **delta_hires:**
%        Neural drive function sampled at 16 Hz (or 'res'), one column per
%        output regressor
%
%   **hrf:**
%        The unit-history (isolated-event) HRF, sampled at 16 Hz
%
%   **out:**
%        Structure with fields:
%          .x_hires       stimulus-history index x(t), one column per
%                         history group (event type)
%          .x_events      cell array, x at each event sample, per regressor
%          .hrf_fcn       function handle: hrf_fcn(x) returns the HRF used
%                         for history index x, sampled at res Hz
%          .hrf_canonical canonical spm_hrf for reference, sampled at res
%          .lambda        decay rate actually used
%          .discount_fcn  handle to w(dt)
%          .res, .TR, .params, .dotimediscount
%
% :Examples:
% ::
%
%    % ------------------------------------------------------------------
%    % Example 1: five convolution models on a simple event sequence
%    % A 4-event 1 Hz burst, an isolated event, then a closely spaced pair.
%    % ------------------------------------------------------------------
%    TR = 1;
%    condf = [1 1 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0]';
%    ons    = {(find(condf) - 1) * TR};
%    runlen = length(condf) * TR;
%
%    Xlin  = conv(condf, spm_hrf(TR) ./ max(spm_hrf(TR)));
%    Xlin  = Xlin(1:length(condf));                       % linear
%    Xmc   = modifiedconv(TR, condf);                     % height only
%    Xvnl  = modifiedconv_wager2005(TR, condf);           % original VNL
%    Xnew0 = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'notimediscount');
%    Xnew  = onsets2fmridesign_nonlin(ons, TR, runlen);   % + time discounting
%
%    t = (0:length(condf)-1)' * TR;
%    colors = seaborn_colors(6);
%    create_figure('Example 1: convolution models'); hold on
%    set(gcf, 'Position', [80 80 900 480]);
%
%    % Event markers along the bottom
%    plot([t(condf > 0) t(condf > 0)]', repmat([-.16; -.04], 1, sum(condf > 0)), ...
%        'Color', [.35 .35 .35], 'LineWidth', 2.5);
%
%    % modifiedconv_wager2005 and the count-history model coincide; the thick
%    % pale band under the thin dashed line shows that.
%    h = [];
%    h(1) = plot(t, Xlin,          'Color', [.55 .55 .55],   'LineWidth', 2.5, 'LineStyle', ':');
%    h(2) = plot(t, Xmc,           'Color', colors{1},       'LineWidth', 2,   'LineStyle', '--');
%    h(3) = plot(t, Xvnl,          'Color', [colors{2} .45], 'LineWidth', 6);
%    h(4) = plot(t, Xnew0(:, 1),   'Color', colors{4},       'LineWidth', 1.5, 'LineStyle', '--');
%    h(5) = plot(t, Xnew(:, 1),    'Color', colors{6},       'LineWidth', 3);
%
%    legend(h, {'Linear convolution' 'modifiedconv (height only)' ...
%        'modifiedconv\_wager2005 (original VNL)' 'nonlin, count history' ...
%        'nonlin, time-discounted'}, 'Location', 'NorthEast'); legend boxoff
%    xlabel('Time (s)'); ylabel('Predicted response (isolated event = 1)');
%    title('Nonlinear saturation: burst, isolated event, pair');
%    set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
%
%    % Note the pair at 20-21 s: the count-based index divides recent history
%    % by the distance back to the unrelated event at 10 s and so predicts no
%    % saturation, while the discounted model correctly gives the second
%    % event of the pair full x = 2 saturation.
%
%    % ------------------------------------------------------------------
%    % Example 2a: HRF shape as a function of the history index x
%    % (as in vnl_equations_plot_2026.m)
%    % ------------------------------------------------------------------
%    TRhi = 0.1;
%    [~, ~, ~, ~, out] = onsets2fmridesign_nonlin({0}, TRhi, 60);
%    xvals = 1:10;
%    colors = seaborn_colors(length(xvals) + 1);
%
%    create_figure('Example 2a: HRF vs history index'); hold on
%    set(gcf, 'Position', [80 80 760 500]);
%    for i = 1:length(xvals)
%        hrf_i = out.hrf_fcn(xvals(i));
%        tt = (0:length(hrf_i) - 1)' ./ out.res;
%        lw = 2; if xvals(i) == 1, lw = 4; end
%        plot(tt, hrf_i, 'Color', colors{i + 1}, 'LineWidth', lw);
%    end
%    set(gca, 'XLim', [0 15], 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
%    xlabel('Time from event (s)'); ylabel('Response (isolated event = 1)');
%    title('HRF vs stimulus history index x (1 Hz train)');
%    legend(arrayfun(@(x) sprintf('x = %d', x), xvals, 'UniformOutput', false), ...
%        'Location', 'NorthEast'); legend boxoff
%
%    % ------------------------------------------------------------------
%    % Example 2b: evoked response to one event, varying time since the
%    % preceding event. This is what time discounting adds.
%    % ------------------------------------------------------------------
%    TRhi = 0.1; lags = [1 2 3 4 5 7];
%    colors = seaborn_colors(length(lags) + 1);
%
%    Xsolo = onsets2fmridesign_nonlin({0}, TRhi, 60);
%    tsolo = (0:size(Xsolo, 1) - 1)' * TRhi;
%
%    create_figure('Example 2b: HRF vs time since last event'); hold on
%    set(gcf, 'Position', [80 80 760 500]);
%    h = plot(tsolo, Xsolo(:, 1), 'k', 'LineWidth', 4);
%    legstr = {'isolated event'};
%
%    for i = 1:length(lags)
%        X = onsets2fmridesign_nonlin({[0; lags(i)]}, TRhi, 60);
%        % History runs only backwards in time, so the first event's response
%        % is unchanged by the second; subtracting it isolates the test event
%        Xtest = X(:, 1) - Xsolo(:, 1);
%        tt = (0:size(X, 1) - 1)' * TRhi - lags(i);   % align on the test event
%        h(end+1) = plot(tt, Xtest, 'Color', colors{i + 1}, 'LineWidth', 2);
%        legstr{end+1} = sprintf('prior event %g s earlier', lags(i));
%    end
%
%    set(gca, 'XLim', [0 20], 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
%    xlabel('Time from event (s)'); ylabel('Response (isolated event = 1)');
%    title('Evoked response vs time since preceding stimulus');
%    legend(h, legstr, 'Location', 'NorthEast'); legend boxoff
%
%    % The discount function itself
%    [~, ~, ~, ~, out] = onsets2fmridesign_nonlin({0}, 0.1, 20);
%    colors = seaborn_colors(3);
%    create_figure('Example 2b: discount function'); hold on
%    set(gcf, 'Position', [80 80 700 460]);
%    dt = linspace(0, 9, 400);
%    plot(dt, out.discount_fcn(dt), 'Color', colors{3}, 'LineWidth', 3);
%    plot([5 5], [0 .1], 'k:', 'LineWidth', 1.5);
%    plot([0 5], [.1 .1], 'k:', 'LineWidth', 1.5);
%    plot([7 7], [0 out.discount_fcn(7)], 'k--', 'LineWidth', 1.5);
%    text(5.2, .14, '10% at 5 s (Miezin et al. 2000)', 'FontSize', 12);
%    text(7.2, .05, 'cutoff', 'FontSize', 12);
%    xlabel('Time since preceding stimulus (s)'); ylabel('Weight w(\Deltat)');
%    title('Exponential discounting of nonlinear saturation');
%    set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
%
%    % ------------------------------------------------------------------
%    % Example 3: dense random event-related design, standard GLM design
%    % matrix construction vs the nonlinear model
%    % ------------------------------------------------------------------
%    TR = 1; runlen = 240;
%    rng(2026);
%    ons = create_random_onsets(runlen, 2, [.25 .25], 1);
%
%    % Impulse events (drop the duration column) so that all five models see
%    % the same input: modifiedconv and modifiedconv_wager2005 take a
%    % TR-sampled indicator vector and cannot represent epochs at all.
%    ons = cellfun(@(x) x(:, 1), ons, 'UniformOutput', false);
%
%    Xlin = onsets2fmridesign(ons, TR, runlen);
%    [Xnl, ~, ~, ~, out] = onsets2fmridesign_nonlin(ons, TR, runlen);
%    Xnl0 = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'notimediscount');
%
%    condf = zeros(runlen / TR, 1);
%    condf(1 + round(ons{1}(:, 1) ./ TR)) = 1;
%    Xmc  = modifiedconv(TR, condf);
%    Xvnl = modifiedconv_wager2005(TR, condf);
%
%    t = (0:size(Xlin, 1) - 1)' * TR;
%    colors = seaborn_colors(6);
%    create_figure('Example 3: dense design', 2, 1);
%    set(gcf, 'Position', [80 80 1000 720]);
%
%    subplot(2, 1, 1); hold on
%    plot(t, Xlin(:, 1), 'Color', [.55 .55 .55],   'LineWidth', 2.5, 'LineStyle', ':');
%    plot(t, Xmc,        'Color', colors{1},       'LineWidth', 1.5, 'LineStyle', '--');
%    plot(t, Xvnl,       'Color', [colors{2} .45], 'LineWidth', 5);
%    plot(t, Xnl0(:, 1), 'Color', colors{4},       'LineWidth', 1.5, 'LineStyle', '--');
%    plot(t, Xnl(:, 1),  'Color', colors{6},       'LineWidth', 2.5);
%    set(gca, 'XLim', [0 120], 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
%    legend({'onsets2fmridesign (linear)' 'modifiedconv' 'modifiedconv\_wager2005' ...
%        'nonlin, count history' 'nonlin, time-discounted'}, 'Location', 'NorthEast');
%    legend boxoff
%    ylabel('Response'); title('Event type 1, first 120 s');
%
%    subplot(2, 1, 2); hold on
%    plot(t, Xlin(:, 1), 'Color', [.55 .55 .55], 'LineWidth', 2.5, 'LineStyle', ':');
%    plot(t, Xlin(:, 2), 'Color', [.75 .75 .75], 'LineWidth', 2.5, 'LineStyle', ':');
%    plot(t, Xnl(:, 1),  'Color', colors{6},     'LineWidth', 2.5);
%    plot(t, Xnl(:, 2),  'Color', colors{3},     'LineWidth', 2.5);
%    set(gca, 'XLim', [0 120], 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
%    legend({'linear, type 1' 'linear, type 2' 'nonlinear, type 1' ...
%        'nonlinear, type 2'}, 'Location', 'NorthEast'); legend boxoff
%    xlabel('Time (s)'); ylabel('Response');
%    title('Both event types (history is within-type only)');
%
%    fprintf('Mean history index x over events: %4.2f (max %4.2f)\n', ...
%        mean(cat(1, out.x_events{:})), max(cat(1, out.x_events{:})));
%    fprintf('Peak amplitude, linear vs nonlinear: %4.2f vs %4.2f\n', ...
%        max(Xlin(:, 1)), max(Xnl(:, 1)));
%
%    % ------------------------------------------------------------------
%    % Other common calls
%    % ------------------------------------------------------------------
%    % Epochs of 3 s, with no post-stimulus undershoot
%    X = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'durs', 3, 'noundershoot');
%
%    % Variance inflation for the nonlinear design
%    create_figure('vifs'); getvif(Xnl, 0, 'plot');
%
%    % One regressor per trial, with history carried across trials
%    X = onsets2fmridesign_nonlin(ons, TR, runlen, [], 'singletrial');
%
% :References:
%   Wager, T. D., Hernandez, L., Vasquez, A., Nichols, T., and Noll, D. C.
%   (2005). Accounting for nonlinear BOLD effects in fMRI: Parameter
%   estimates and model for accurate prediction in variable-duration blocked
%   and rapid event-related studies. Neuroimage, 25(1), 206-218.
%
%   Miezin, F. M., Maccotta, L., Ollinger, J. M., Petersen, S. E., and
%   Buckner, R. L. (2000). Characterizing the hemodynamic response: effects
%   of presentation rate, sampling procedure, and the possibility of
%   ordering brain activity based on relative timing. NeuroImage, 11, 735-759.
%
% :See also: onsets2fmridesign, modifiedconv_wager2005, modifiedconv,
%            plotDesign, getvif

% ..
%    Tor Wager and CANlab, 2026
%
%    Programmers' notes:
%    This function was originally split into a convolution engine
%    (conv_vnl_nonlin.m) and this wrapper. The engine had exactly one
%    caller, and each half carried its own copy of the onset-normalization
%    and duration-appending helpers, so the two were merged here in 2026.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Positional optional inputs are [len] and [HRF]. If the first optional
% argument is a string, the caller passed keywords only.
if ~isempty(varargin) && ischar(varargin{1})
    n_fixed_args = 0;
else
    n_fixed_args = 2;
end

len = [];
if n_fixed_args > 0 && ~isempty(varargin)
    len = varargin{1};

    if ~isempty(len) && ~(isnumeric(len) && isscalar(len) && len > 0)
        error('Run length (3rd argument) must be a positive scalar in seconds, or [].');
    end
end

if n_fixed_args > 0 && length(varargin) > 1 && ~isempty(varargin{2})
    error(['onsets2fmridesign_nonlin does not accept a custom HRF or basis set. ' ...
        'The nonlinear model is parameterized through spm_hrf. Pass [] instead.']);
end

varargin(1:min(n_fixed_args, length(varargin))) = [];

% Parse special command keywords and remove them before inputParser

doplot = false;
wh = strcmpi(varargin, 'plot') | strcmpi(varargin, 'doplot');
if any(wh), doplot = true; varargin(wh) = []; end

doverbose = false;
wh = strcmpi(varargin, 'verbose');
if any(wh), doverbose = true; varargin(wh) = []; end
wh = strcmpi(varargin, 'noverbose');
if any(wh), doverbose = false; varargin(wh) = []; end

donoundershoot = false;
wh = strcmpi(varargin, 'noundershoot');
if any(wh), donoundershoot = true; varargin(wh) = []; end

dotimediscount = true;
wh = strcmpi(varargin, 'notimediscount') | strcmpi(varargin, 'nodiscount') | strcmpi(varargin, 'original');
if any(wh), dotimediscount = false; varargin(wh) = []; end

donorm = true;
wh = strcmpi(varargin, 'nonorm');
if any(wh), donorm = false; varargin(wh) = []; end

dosingletrial = false;
wh = strcmpi(varargin, 'singletrial');
if any(wh), dosingletrial = true; varargin(wh) = []; end

input_is_condf = false;
wh = strcmpi(varargin, 'condf');
if any(wh), input_is_condf = true; varargin(wh) = []; end

input_is_delta_hires = false;
wh = strcmpi(varargin, 'deltahires') | strcmpi(varargin, 'delta_hires');
if any(wh), input_is_delta_hires = true; varargin(wh) = []; end

% Options that belong to onsets2fmridesign but have no meaning here
wh = strcmpi(varargin, 'nonlinsaturation') | strcmpi(varargin, 'nononlin');
if any(wh), varargin(wh) = []; end          % nonlinearity is intrinsic

wh = strcmpi(varargin, 'norm') | strcmpi(varargin, 'noampscale');
if any(wh)
    warning('''%s'' applies to basis sets and is ignored by onsets2fmridesign_nonlin.', varargin{find(wh, 1)});
    varargin(wh) = [];
end

wh = strcmpi(varargin, 'parametric_standard') | strcmpi(varargin, 'parametric_singleregressor');
if any(wh)
    error(['Parametric modulators are not implemented in onsets2fmridesign_nonlin. ' ...
        'Use onsets2fmridesign.']);
end

if any(strcmpi(varargin, 'customneural'))
    error(['Custom neural response functions are not implemented as a keyword. ' ...
        'Build the drive function yourself and pass it with the ''deltahires'' keyword, ' ...
        'or use onsets2fmridesign.']);
end

% Use inputParser to parse key/value pairs

valfcn_eq = @(x) isempty(x) || isa(x, 'function_handle');
valfcn_pos = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', 'positive'});

p = inputParser;
p.addRequired('ons', @(x) isnumeric(x) || iscell(x));
p.addRequired('TR', valfcn_pos);

p.addParameter('durs', [], @(x) isempty(x) || isnumeric(x) || iscell(x));
p.addParameter('res', 16, valfcn_pos);
p.addParameter('decay5', 0.10, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<', 1}));
p.addParameter('lambda', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x > 0));
p.addParameter('maxhistory', 7, valfcn_pos);
p.addParameter('crossweights', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('hrfparams', [6 16 1 1 6 0 32], @(x) isnumeric(x) && numel(x) == 7);
p.addParameter('xstep', 0.02, valfcn_pos);   % quantization of x for HRF caching
p.addParameter('maxx', 20, valfcn_pos);      % ceiling on the history index

% Wager et al. 2005 equations. x = 1 is an isolated event.
p.addParameter('heighteq',  @(x)   1.7141 .* exp(-2.1038 .* x) + 0.4932 .* exp(-0.0770 .* x), valfcn_eq);
p.addParameter('delayeq',   @(x) -13.4097 .* exp(-1.0746 .* x) + 4.8733 .* exp(-0.1979 .* x), valfcn_eq);
p.addParameter('ttopeakeq', @(x)  37.5445 .* exp(-2.6760 .* x) - 3.2046 .* exp(-0.2120 .* x) + 5.6344, valfcn_eq);
p.addParameter('uonseteq',  [], valfcn_eq);

p.parse(ons, TR, varargin{:});

ARGS = p.Results;
res         = ARGS.res;
durs        = ARGS.durs;
maxhistory  = ARGS.maxhistory;
hrfparams   = ARGS.hrfparams;
xstep       = ARGS.xstep;
maxx        = ARGS.maxx;
heighteq    = ARGS.heighteq;
delayeq     = ARGS.delayeq;
ttopeakeq   = ARGS.ttopeakeq;
uonseteq    = ARGS.uonseteq;

if donoundershoot, hrfparams(5) = Inf; end

if input_is_delta_hires && dosingletrial
    error('''singletrial'' needs an onset list; it cannot be combined with ''deltahires''.');
end

% Exponential discount rate. w(1) = 1 by construction; lambda is chosen so
% that w(5) = decay5.
if isempty(ARGS.lambda)
    lambda = -log(ARGS.decay5) ./ 4;
else
    lambda = ARGS.lambda;
end

discount_fcn = @(dt) exp(-lambda .* max(0, dt - 1));

% The HRF actually used at a given stimulus-history index
hrf_scale = 1;
if donorm && ~isempty(heighteq), hrf_scale = 1 ./ heighteq(1); end

hrf_fcn = @(x) hrf_scale .* local_vnl_hrf(1 ./ res, hrfparams, x, ...
    heighteq, delayeq, ttopeakeq, uonseteq);

% -------------------------------------------------------------------------
% Build the high-resolution neural drive function
%
% delta_hist  [samples x n history groups]  drives that generate history
% delta_hires [samples x n regressors]      drives that generate regressors
% x_source    [1 x n regressors]            which history group each uses
%
% These differ only for single-trial models, where each trial gets its own
% regressor but history still accumulates across all trials of the original
% event type.
% -------------------------------------------------------------------------

if input_is_delta_hires

    delta_hist = double(ons);
    if isempty(len), len = size(delta_hist, 1) ./ res; end

    delta_hires = delta_hist;
    x_source = 1:size(delta_hist, 2);
    ons = {};

else

    if input_is_condf
        % Indicator vector sampled at TR, with integer event codes
        condf = double(ons(:));
        if isempty(len), len = length(condf) .* TR; end

        n_conditions = max([0; condf]);
        onscell = cell(1, n_conditions);
        for i = 1:n_conditions
            onscell{i} = (find(condf == i) - 1) .* TR;   % element 1 is time 0
        end
        ons = onscell;
    end

    ons = local_check_onsets(ons);

    if ~isempty(durs)
        ons = local_add_durs(ons, durs);
    end

    % Length of the run
    if ~isempty(len)
        for i = 1:length(ons)
            past_end = ons{i}(:, 1) > len;
            if any(past_end)
                fprintf(['Warning! onsets2fmridesign_nonlin found %d onsets after the ' ...
                    'end of the session and is omitting them.\n'], sum(past_end));
            end
            ons{i}(past_end, :) = [];
        end
        ons = local_check_onsets(ons);

        len_hires = ceil(len .* res);
    else
        % Grow to accommodate onsets; round up to a whole number of TRs
        maxons = max(cat(1, ons{:}));
        len_hires = round(maxons(1) .* res);
        len_hires = ceil(ceil(len_hires ./ (res .* TR)) .* res .* TR);
        len = len_hires ./ res;
    end

    % History groups: always the original event types
    delta_hist = local_build_drive(ons, len_hires, res);

    if dosingletrial
        % One regressor per trial, but keep the history grouping above
        sons = {};
        x_source = [];

        for i = 1:length(ons)
            for j = 1:size(ons{i}, 1)
                sons{end+1} = ons{i}(j, :);     %#ok<AGROW>
                x_source(end+1) = i;            %#ok<AGROW>
            end
        end

        ons = sons;
        delta_hires = local_build_drive(ons, len_hires, res);

    else
        delta_hires = delta_hist;
        x_source = 1:length(ons);
    end

end

len_hires = size(delta_hires, 1);
n_regressors = size(delta_hires, 2);
n_groups = size(delta_hist, 2);

% -------------------------------------------------------------------------
% Stimulus-history index x(t) for each event type
% -------------------------------------------------------------------------

crossweights = ARGS.crossweights;
if isempty(crossweights)
    crossweights = eye(n_groups);    % no cross-event-type influence
elseif ~isequal(size(crossweights), [n_groups n_groups])
    error('crossweights must be %d x %d (one row/column per event type).', n_groups, n_groups);
end

if dotimediscount
    % Discounted history: convolve each event type's drive with the decay
    % kernel, offset by one sample so a stimulus never saturates itself.
    nback = max(1, round(maxhistory .* res));
    kernel = [0; discount_fcn((1:nback)' ./ res)];

    hsum = zeros(len_hires, n_groups);
    for i = 1:n_groups
        tmp = conv(delta_hist(:, i), kernel);
        hsum(:, i) = tmp(1:len_hires);
    end

    x_hires = 1 + hsum * crossweights';

else
    % Original count-based index from modifiedconv_wager2005: number of
    % same-type events in the last 12 s, scaled by their density in
    % TR-sized bins.
    if any(delta_hist(delta_hist ~= 0) ~= 1)
        warning(['The original count-based history index counts nonzero neural samples, ' ...
            'so it is only meaningful for impulse events. You have entered epochs or ' ...
            'graded amplitudes; every high-resolution sample of an epoch will be counted ' ...
            'as a separate event. Use time discounting (the default) instead.']);
    end

    x_hires = local_count_history(delta_hist, TR, res);

    if ~isequal(crossweights, eye(n_groups))
        warning('crossweights are ignored when time discounting is off.');
    end
end

x_hires = max(1, min(maxx, x_hires));

% -------------------------------------------------------------------------
% Convolve each event with its own history-dependent HRF
% -------------------------------------------------------------------------

hrf_len = round(hrfparams(7) .* res) + 1;

% HRF cache, keyed on the quantized history index
hrf_cache = containers.Map('KeyType', 'double', 'ValueType', 'any');

model_hires = zeros(len_hires + hrf_len, n_regressors);
x_events = cell(1, n_regressors);

for i = 1:n_regressors

    s = delta_hires(:, i);
    whs = find(s ~= 0);

    xi = x_hires(whs, x_source(i));
    x_events{i} = xi;

    xkey = round((xi - 1) ./ xstep);

    for a = 1:length(whs)

        key = xkey(a);

        if hrf_cache.isKey(key)
            h = hrf_cache(key);
        else
            h = hrf_fcn(1 + key .* xstep);
            h = [h; zeros(max(0, hrf_len - length(h)), 1)];   % pad to fixed length
            h = h(1:hrf_len);
            hrf_cache(key) = h;
        end

        t0 = whs(a);
        model_hires(t0:t0 + hrf_len - 1, i) = model_hires(t0:t0 + hrf_len - 1, i) + s(t0) .* h;

    end

end

model_hires = model_hires(1:len_hires, :);

% -------------------------------------------------------------------------
% Downsample to TR
% (matches getPredictors.m, so output is aligned with onsets2fmridesign)
% -------------------------------------------------------------------------

n = size(model_hires, 1);
tvec = (1:n)';

nt = round(len ./ TR);
ds_rate = n ./ nt;
xq = 1:ds_rate:n;

X = zeros(length(xq), n_regressors);
for i = 1:n_regressors
    X(:, i) = interp1(tvec, model_hires(:, i), xq, 'linear', 'extrap');
end

if size(X, 1) > nt, X = X(1:nt, :); end

nvols = size(X, 1);

X(:, end+1) = 1;    % intercept, as in onsets2fmridesign

% -------------------------------------------------------------------------
% TR-resolution indicator matrix (for HRF shape estimation), as in
% onsets2fmridesign
% -------------------------------------------------------------------------

delta = cell(1, length(ons));

for i = 1:length(ons)

    delta{i} = zeros(nvols, 1);

    tmp = 1 + round(ons{i}(:, 1) ./ TR);     % time 0 is first element
    tmp(tmp > nvols) = [];

    delta{i}(tmp) = 1;

    repeats = tmp(diff([tmp; Inf]) == 0);
    delta{i}(repeats) = delta{i}(repeats) + 1;

end

% -------------------------------------------------------------------------
% Remaining outputs, verbose report, and optional plot
% -------------------------------------------------------------------------

hrf = hrf_fcn(1);

out = struct();
out.x_hires = x_hires;
out.x_events = x_events;
out.hrf_fcn = hrf_fcn;
out.hrf_canonical = spm_hrf(1 ./ res, hrfparams);
out.hrf_canonical = out.hrf_canonical ./ max(out.hrf_canonical);
out.lambda = lambda;
out.discount_fcn = discount_fcn;
out.res = res;
out.TR = TR;
out.params = hrfparams;
out.dotimediscount = dotimediscount;

if doverbose
    fprintf('onsets2fmridesign_nonlin: %d regressor(s), %3.1f s at TR = %3.2f s (%d samples)\n', ...
        n_regressors, len, TR, nvols);
    if dotimediscount
        fprintf('  Time-discounted history: lambda = %3.4f, w(2 s) = %3.3f, w(5 s) = %3.3f, cutoff %3.1f s\n', ...
            lambda, discount_fcn(2), discount_fcn(5), maxhistory);
    else
        fprintf('  Original count-based history index (no time discounting)\n');
    end
    fprintf('  History index x: max = %3.2f, mean over events = %3.2f\n', ...
        max([1; cat(1, x_events{:})]), mean([1; cat(1, x_events{:})]));
end

if doplot
    create_figure('onsets2fmridesign_nonlin');
    t = (0:nvols - 1)' .* TR;
    colors = seaborn_colors(max(3, n_regressors));
    hold on
    for i = 1:n_regressors
        plot(t, X(:, i), 'Color', colors{i}, 'LineWidth', 2);
    end
    xlabel('Time (s)'); ylabel('Predicted response');
    set(gca, 'FontSize', 14, 'Box', 'off', 'TickDir', 'out');
end

end % main function



% =========================================================================
% Sub-functions
% =========================================================================

function drive = local_build_drive(ons, len_hires, res)
% High-resolution neural drive: unit amplitude per impulse event, or 1/res
% per sample of an epoch, so a 1 s epoch delivers one unit of drive. This
% matches onsets2fmridesign's neural function exactly.

drive = zeros(len_hires, length(ons));

for i = 1:length(ons)

    has_durations = size(ons{i}, 2) > 1;

    for j = 1:size(ons{i}, 1)

        fs = round(ons{i}(j, 1) .* res) + 1;    % time 0 is element 1
        if fs > len_hires, continue, end

        if ~has_durations || ons{i}(j, 2) <= 0
            drive(fs, i) = drive(fs, i) + 1;
        else
            es = min(len_hires, round(fs + ons{i}(j, 2) .* res));
            drive(fs:es, i) = 1 ./ res;
        end

    end

end

end



function h = local_vnl_hrf(dt, hrfparams, x, heighteq, delayeq, ttopeakeq, uonseteq)
% HRF for stimulus-history index x, sampled every dt seconds, scaled so that
% its peak equals heighteq(x).

trialp = hrfparams;

if ~isempty(ttopeakeq), trialp(1) = ttopeakeq(x); end
if ~isempty(uonseteq),  trialp(2) = uonseteq(x);  end
if ~isempty(delayeq),   trialp(6) = delayeq(x);   end

h = spm_hrf(dt, trialp);

mx = max(h);
if mx <= 0
    % Degenerate parameters; return a flat response rather than dividing by zero
    h = zeros(size(h));
    return
end

h = h ./ mx;

if ~isempty(heighteq), h = heighteq(x) .* h; end

end



function x_hires = local_count_history(delta_hist, TR, res)
% Original modifiedconv_wager2005 history index, evaluated at the
% high-resolution event samples: number of same-type events within the
% preceding 12 s (including the current one), scaled by their density in
% TR-sized bins, floored at 1.

[len_hires, k] = size(delta_hist);
x_hires = ones(len_hires, k);

nback = max(1, round(12 .* res));      % 12 s window
binsize = max(1, res .* TR);           % samples per TR

for i = 1:k

    whs = find(delta_hist(:, i) ~= 0);
    lo = 1;

    for a = 1:length(whs)

        t = whs(a);
        while whs(lo) < t - nback, lo = lo + 1; end

        stimpos = a - lo + 1;                                  % events in window, incl. current
        nbins = round((t - whs(lo)) ./ binsize) + 1;           % TR-sized bins spanned
        dens = stimpos ./ nbins;

        x_hires(t, i) = max(1, stimpos .* dens);

    end

end

end



function ons = local_check_onsets(ons)
% Enforce a cell array of [n_events x 1] or [n_events x 2] onset matrices

if ~iscell(ons), ons = {ons}; end

for i = 1:length(ons)

    if isempty(ons{i})
        error('Onsets for event type %d are empty. Check inputs.', i);
    end

    % A row vector of onsets with no durations is a common input; transpose.
    if size(ons{i}, 1) == 1 && size(ons{i}, 2) ~= 2
        ons{i} = ons{i}';
    end

    if any(ons{i}(:, 1) < 0)
        error('Illegal onsets (negative values) for event type %d.', i);
    end

end

end



function ons = local_add_durs(ons, durs)
% Append durations (in SECONDS) as column 2 of each onset matrix, matching
% onsets2fmridesign.

if all(cellfun(@(x) size(x, 2) > 1, ons))
    warning('Onsets input already contains durations. ''durs'' input will be ignored.');
    return
end

if iscell(durs)

    if length(durs) ~= length(ons)
        error('''durs'' cell array has %d cells but there are %d event types.', length(durs), length(ons));
    end

    for i = 1:length(ons)
        if numel(durs{i}) ~= size(ons{i}, 1)
            error('''durs''{%d} has %d durations but event type %d has %d events.', ...
                i, numel(durs{i}), i, size(ons{i}, 1));
        end
        ons{i} = [ons{i} durs{i}(:)];
    end

elseif isscalar(durs)
    for i = 1:length(ons), ons{i} = [ons{i} repmat(durs, size(ons{i}, 1), 1)]; end

elseif length(durs) == length(ons)
    for i = 1:length(ons), ons{i} = [ons{i} repmat(durs(i), size(ons{i}, 1), 1)]; end

else
    warning('''durs'' input format/length is unrecognized. It will not be used.');

end

end
