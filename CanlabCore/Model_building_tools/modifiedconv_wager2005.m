function model = modifiedconv_wager2005(tr, indicator_vec, varargin)
% Modified convolution adjusting response height, onset delay, and time to
% peak as a function of stimulus history (Wager et al. 2005 equations)
%
% :Usage:
% ::
%
%     model = modifiedconv_wager2005(tr, indicator_vec, ['height', heighteq], ['delay', delayeq], ['ttopeak', peakeq], ['uonset', uonseteq])
%
% This is the original ("VNL") implementation, retained for reference and
% for regression-testing newer code. It works on an indicator vector
% sampled at the TR, treats stimulus history as a *count* of recent events
% (scaled by a crude density correction), and is only strictly valid for
% 1 s interstimulus intervals. For a version that handles onsets/durations
% in seconds, multiple event types, and time-based discounting of
% nonlinearity, see onsets2fmridesign_nonlin.m.
%
% :Inputs:
%
%   **tr:**
%        repetition time (sampling rate) of scanning, in seconds
%
%   **indicator_vec:**
%        condition function: a vector of zeros and integers, where integer
%        i marks an onset of event type i. One output column per event type.
%
% :Optional Inputs:
%
%   **'height':**
%        followed by a function handle mapping stimulus-history index x to
%        response height. Enter [] to disable height modulation.
%
%   **'delay':**
%        followed by a function handle mapping x to onset delay in sec
%        (spm_hrf parameter 6). Enter [] to disable.
%
%   **'ttopeak':**
%        followed by a function handle mapping x to time to peak in sec
%        (spm_hrf parameter 1). Enter [] to disable.
%
%   **'uonset':**
%        followed by a function handle mapping x to undershoot onset delay
%        in sec (spm_hrf parameter 2). Off by default.
%
%   **'verbose':**
%        print the equations in use
%
% :Outputs:
%
%   **model:**
%        [length(indicator_vec) x n event types] predictors, sampled at tr,
%        normalized so that an isolated (unit-history) event has peak 1.
%
% :Examples:
% ::
%
%    indicator_vec = [1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]';
%    X = modifiedconv_wager2005(2, indicator_vec);
%    figure; plot(X)
%    X2 = conv(indicator_vec, spm_hrf(2)./max(spm_hrf(2)));
%    hold on; plot(X2(1:length(X)), 'r'); legend({'Modified' 'Linear'})
%
% :References:
%   Wager, T. D., Hernandez, L., Vasquez, A., Nichols, T., and Noll, D.
%   C. (2005). Accounting for nonlinear BOLD effects in fMRI: Parameter
%   estimates and model for accurate prediction in variable-duration blocked
%   and rapid event-related studies. Neuroimage.
%
% :See also: modifiedconv, onsets2fmridesign_nonlin, onsets2fmridesign
%
% ..
%    06/20/01 Tor Wager
%    2026     Maintenance pass (see Programmers' notes at end of file)
% ..

% ---------------------------------------------------------------------
% * defaults
% ---------------------------------------------------------------------
height = 1; delay = 0; peak = 6; uonset = 16;
dispers = 1; udisp = 1; rtou = 6; klength = 32;

p = [peak uonset dispers udisp rtou delay klength];

% These equations use parameters estimated from Wager et al. 2005.
% x is the stimulus-history index: x = 1 is an isolated event, x = 2 is the
% second event of a 1 Hz train, and so on.
heighteq = @(x)   1.7141 .* exp(-2.1038 .* x) + 0.4932 .* exp(-0.0770 .* x);
delayeq  = @(x) -13.4097 .* exp(-1.0746 .* x) + 4.8733 .* exp(-0.1979 .* x);
peakeq   = @(x)  37.5445 .* exp(-2.6760 .* x) - 3.2046 .* exp(-0.2120 .* x) + 5.6344;
uonseteq = [];

doverbose = false;

% ---------------------------------------------------------------------
% * set up arguments
% ---------------------------------------------------------------------

ind = 1;
while ind <= length(varargin)

    if ~ischar(varargin{ind})
        ind = ind + 1;
        continue
    end

    switch varargin{ind}

        case 'verbose'
            doverbose = true;
            ind = ind + 1;
            continue

        case {'height', 'delay', 'ttopeak', 'uonset'}
            if ind + 1 > length(varargin)
                error('Option ''%s'' must be followed by a function handle or [].', varargin{ind});
            end

        otherwise
            error('Unknown input string option: %s', varargin{ind});
    end

    switch varargin{ind}
        case 'height',  heighteq = varargin{ind+1};
        case 'delay',   delayeq  = varargin{ind+1};
        case 'ttopeak', peakeq   = varargin{ind+1};
        case 'uonset',  uonseteq = varargin{ind+1};
    end

    ind = ind + 2;
end

if doverbose
    fprintf('modifiedconv_wager2005: height   = %s\n', eqstring(heighteq));
    fprintf('modifiedconv_wager2005: delay    = %s\n', eqstring(delayeq));
    fprintf('modifiedconv_wager2005: ttopeak  = %s\n', eqstring(peakeq));
    fprintf('modifiedconv_wager2005: uonset   = %s\n', eqstring(uonseteq));
end

% ---------------------------------------------------------------------
% * build predictors
% ---------------------------------------------------------------------

indicator_vec = double(indicator_vec(:));

myzeros = zeros(length(indicator_vec), 1);
mylen = length(myzeros);

% Number of elements of history to consider (12 sec). Must be an integer
% number of samples, or the indexing below is not a valid subscript.
numels = max(1, round(12 ./ tr));

n_conditions = max([0; indicator_vec]);
model = zeros(mylen, n_conditions);

for i = 1:n_conditions
    % For each integer in indicator_vec

    reg = myzeros;

    for j = 1:length(indicator_vec)
        % for each time point

        if indicator_vec(j) == i

            trialdelta = myzeros;
            trialdelta(j) = 1;
            trialp = p;

            % figure out how many of the same type came before
            % 1 is "first stim in sequence"

            myc = indicator_vec(j-min(numels, j)+1 : j);

            stimpos = sum(myc == i);
            tmp = find(myc == i);

            % Empty-to-filled ratio within the window, counting from the
            % first event of THIS type. This is a scaling factor for
            % stimulus position based on history density: if prior history
            % is not full of stimulation, the 'position' estimate is
            % shifted toward 1, but cannot go below 1.
            dens = stimpos ./ length(myc(tmp(1):end));

            x = max(1, stimpos .* dens);

            if ~isempty(heighteq)
                height = heighteq(x);
            end

            if ~isempty(delayeq)
                trialp(6) = delayeq(x);
            end

            if ~isempty(peakeq)
                trialp(1) = peakeq(x);
            end

            if ~isempty(uonseteq)
                trialp(2) = uonseteq(x);
            end

            trialhrf = spm_hrf(tr, trialp);
            trialhrf = height .* (trialhrf ./ max(trialhrf));

            mytrialpred = conv(trialhrf, trialdelta);
            mytrialpred = mytrialpred(1:mylen);

            reg = reg + mytrialpred;

        end  % time point
    end % integer

    % normalize so that height of unit response (stimpos == 1) is 1
    if ~isempty(heighteq)
        model(:, i) = reg ./ heighteq(1);
    else
        model(:, i) = reg;
    end

end

end % main function



function s = eqstring(eq)

if isempty(eq)
    s = '(none)';
elseif isa(eq, 'function_handle')
    s = func2str(eq);
else
    s = class(eq);
end

end

% ..
%    Programmers' notes:
%
%    2026 maintenance pass. The original code already applied all three
%    fitted equations (height -> amplitude, delayeq -> spm_hrf p(6) onset
%    delay, peakeq -> spm_hrf p(1) time to peak); that behavior is
%    unchanged for single-event-type, integer-TR input. Changes made:
%
%    1. numels = 12./tr is now rounded. With a fractional TR (e.g. 1.3) the
%       unrounded value was used as a colon-operator subscript, which threw
%       "Integer operands are required for colon operator when used as
%       index" and silently truncated the history window.
%    2. The local variable 'disp' (dispersion) was renamed 'dispers'. It
%       shadowed the built-in disp() for the rest of the function's scope.
%    3. Stimulus-history density is now computed from events of the SAME
%       event type (find(myc == i)) rather than from any nonzero entry.
%       With more than one event type the old version diluted the density
%       correction with other conditions' events, so history for type A
%       depended on type B. No change for single-event-type input.
%    4. 'uonset' now receives the same density-corrected history index as
%       the other three equations, instead of the raw event count.
%    5. inline() (deprecated) replaced with anonymous function handles.
%    6. Option parsing bounds-checks its arguments, errors on unrecognized
%       strings, and no longer echoes the equations to the command window
%       on every call (use 'verbose' for that).
%    7. 'model' is preallocated, so an all-zero indicator vector returns an
%       empty matrix instead of erroring on an undefined output.
%    8. Removed an unused precomputed canonical hrf (dead code).
% ..
