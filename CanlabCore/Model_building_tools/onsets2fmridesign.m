function [X, delta, delta_hires, hrf] = onsets2fmridesign(ons, TR, varargin)
% Builds design matrix X and delta function, given cell array of onset times in s
% One cell per condition per session, e.g., ons{1} = [24 27 29 44]';
%
% :Usage:
% ::
%
%     [X, delta, delta_hires, hrf] = onsets2fmridesign(onsets, TR, [len], [custom hrf or basis set name],[other optional args, i.e., norm])
%
% Summary:
%   - handles multiple conditions
%   - handles custom HRFs and multiple basis functions
%   - handles input event durations, variable-duration stimuli
%   - handles two kinds of parametric modulators
%   - limited custom neural response function
%   - handles nonlinear saturation (see hrf_saturation.m)
%   - Can build single-trial model
%   - not yet: variable-duration parametric
%     modulators
%   - See the code comments for a discussion of absolute scaling of
%     regressors and efficiency.
%
% :Inputs:
%
%   **onsets:**
%        - 1st column is onsets for events in seconds,
%        - 2nd column is optional durations for each event
%        - Enter single condition or cell vector with cells for each condition (each event type).
%
%   **TR:**
%        Repetition time (RT) in seconds
%
% :Optional Inputs: First two are fixed, then keywords:
%
%   **len:**
%        optional length in sec for model, or [] to skip if using additional. 
%        "len" is usually the number of images multiplied by TR, n_images * TR
%
%   **HRF name:**
%        string or actual values
%             1) a string used by spm_get_bf.m
%             2) a custom HRF, sampled in seconds
%             3) or [] to skip
%
%   **'norm':**
%        mean-center, orthogonalize, and L2-norm basis set
%
%   **'parametric_singleregressor':**
%        Parametrically modulate onsets by
%               modulator values, using single regressor with modulated amplitude
%               Enter a cell array with modulator values for each event
%               type, with a column vector (empty cell for no modulation)
%
%   **'parametric_standard':**
%        Parametrically modulate onsets by
%        modulator values, using two regressors per event type - One
%        to model the average response, and one for the
%        mean-centered modulator values
%        of modulator values in each cell
%
%   **'noampscale':**
%        Do not scale HRF to max amplitude = 1; default with SPM
%        basis sets is to scale.  Custom HRF entries are not scaled.
%
%   **'noundershoot':**
%        Do not model undershoot - ONLY when using the canonical HRF
%
%   **'customneural':**
%        Followed by a vector of custom neural values (instead of
%        standard event/epochs), sampled in sec.
%
% Limitation: can only handle two events max within the same TR
%
% :Outputs:
%
%   **X:**
%        model, sampled at TR
%
%   **delta:**
%        indicator mtx, sampled at TR
%
%   **delta_hires:**
%        indicator sampled at 16 Hz
%
%   **hrf:**
%        hemodynamic response function, sampled at 16 Hz
%
% :Examples:
% ::
%
%    X = onsets2fmridesign(ons, TR, size(imgs, 1) .* TR, 'hrf (with time derivative)');
%
%    X = onsets2fmridesign({[0 30 60]' [15 45 90]'}, 1.5, 180, spm_hrf(1), 'parametric_standard', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
%    X = onsets2fmridesign({[0 30 60]' [15 45 90]'}, 1, 180, spm_hrf(1), 'parametric_standard', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
%    X = onsets2fmridesign({[0 30 60]' [15 45 90]'}, 2, 180, spm_hrf(1), 'parametric_standard', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
%    X = onsets2fmridesign({[0 30 60]' [15 45 90]'}, 2, 180, spm_hrf(1), 'parametric_singleregressor', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
%    X = onsets2fmridesign({[0 30 60]' [15 45 90]'}, 1, 180, spm_hrf(1), 'parametric_singleregressor', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
%
%    % Here, spm_hrf(1) is for the canonical HRF in spm.
%
%    X = onsets2fmridesign(ons, 2, size(dat.dat, 2) .* 2, [], [], 'noampscale');
%
%    X2 = onsets2fmridesign(onsp, 2, length(dat.removed_images) .* 2, [], [], 'customneural', report_predictor);
%    plotDesign(ons,[], TR, 'yoffset', 1);
%
% A note on absolute scaling and efficiency:
% Scaling of the response influences efficiency. It does not affect model
% fits or power when the scaling is equated (constant across events in a
% design), but it does affect efficiency simulations.
% 
% Event duration:
%    assumes that the max neural sampling rate is about 1 Hz,
%    which produces responses that do not exceed about 5x the
%    unit, single-event response.  This is a reasonable
%    assumption, and affects only the absolute scaling across
%    ISIs/block lengths. 
% 
% An event duration of about 1 sec produces a response of unit amplitude.
% The sampling resolution is 0.1 sec, so that is the lowest you can go - and it will produce lower response amplitude/less efficient designs, as less neural activity is being delivered.
% If event durations are not entered, then events will have unit amplitude by default.
%
% :See Also: tor_make_deconv_mtx3, [or 2], plotDesign

%
% ..
%    Programmers' notes:
%    Tor Wager, 2 / 24 / 04
%    Modified/updated 2/08 by Tor to add capacity for durations
%    Modified 10/1/11 by Tor to fix custom HRF sampling
%    Modified 3/14/12 by Tor to check variable duration and add noampscale option
%    3/16/12: Tor modified custom HRF to divide by sum, as SPM basis functions
%    are scaled by default.
%    9/3/12: Wani modified len and getPredictors.m to fix a bug when you're
%    using TR with some decimal points (e.g., TR = 1.3)
%    Modified 8/2015 by Tor to clean up documentation and change scaling of
%    regressors overall in better range
%    1/2020: Tor modified to increase tolerance/flexibilty for fractional
%    TRs, add better error reporting for length mismatch
% ..

% ..
%    Defaults
% ..

res = 16;   % resolution, in samples per second

% Other optional inputs after fixed ones
% - - - - - - - - - - - - - - - - - - - -
donorm = 0;
doampscale = 1;
hrfparams = [6 16 1 1 6 0 32];
customneural = []; % custom neural response; neither impulse nor epoch
dosingletrial = 0;
docheckorientation = 1;
dononlinsaturation = 0;

for i = 3:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % reserved keywords
            case 'norm', donorm = 1;
                
            case {'parametric_standard', 'parametric_singleregressor'}
                % pm_mods = varargin{i + 1}; unnecessary, as passed
                % into getPredictors
                
            case 'noampscale', doampscale = 0;
                
            case 'noundershoot', hrfparams(5) = Inf;
                
            case 'customneural', customneural = varargin{i + 1};
                
            case 'singletrial', dosingletrial = 1; docheckorientation = 0;
                
            case 'nonlinsaturation', dononlinsaturation = 1;
                
            case 'nononlin',  dononlinsaturation = false;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Check that ons is a cell array, orientation, and check numbers of events
% -------------------------------------------------------------------------
[ons, n_conditions, n_events, ons_includes_durations] = check_onsets(ons, docheckorientation);
%
% [ons, n_conditions, n_events, ons_includes_durations] = check_onsets(ons, docheckorientation)
% n_conditions = number of event types, length(ons)
% n_events = vector of number of events in each condition

% ----------------------------------------------
% Set up single trial: optional
% ----------------------------------------------

if dosingletrial
    % single trial: break events into separate regressors/onset vectors
    for i = 1:length(ons)
        
        sonsets{i} = cell(1, n_events(i)); % event-specific onsets - init
        
        for j = 1:length(sonsets{i})
            
            sonsets{i}{j} = ons{i}(j, 1);
            if size(ons{i}, 2) > 1   % if durations...
                sonsets{i}{j} = [sonsets{i}{j} ons{i}(j, 2)]; % append duration
            end
            
        end
        
    end % original onset conditions
    
    ons = cat(2, sonsets{:}); % flatten to simply model trials within event type.
    n_conditions = sum(n_events);
    n_events = ones(n_conditions, 1);

end  % end single trial


% Optional: pre-specified length for run, or []
% - - - - - - - - - - - - - - - - - - - -

if ~isempty(varargin) && ~isempty(varargin{1}) % pre-specified length
    % Entering session len, and using len_original here, will downsample using pre-specified length
    
    len_original = varargin{1}; % this line is added by Wani to make this work for TR = 1.3
    len = varargin{1}; % this line is modified by Wani from the next line
    % len = ceil(varargin{1});
    for i = 1:length(ons)

        % make sure nothing goes past end of session
        past_end = round(ons{i}(:,1)) > len;
        if any(past_end)
            fprintf('Warning! onsets2fmridesign found %d onsets after the end of the session and is omitting them.\n', sum(past_end));
        end
        ons{i}(past_end,:) = [];
    end
    len = ceil(len .* res);
    
else
    % No length entered; length is allowed to grow to accomodate onsets
    
    len = round(max(cat(1, ons{:})) * res);
    len = ceil(ceil(len/(res*TR)) * res * TR); % modified by Wani to make this work for TR = 1.3
    len = len(1);
    
    len_original = []; % Use 'dsrate' to downsample
end

% len was in sec, now in sec*res
cf = zeros(len, 1);

% Optional: basis set name, or []
% - - - - - - - - - - - - - - - - - - - -
if length(varargin) > 1 && ~isempty(varargin{2})
    if ischar(varargin{2})
        % basis set: Any of the named basis sets from SPM  spm_get_bf
        
        % These need 'order' parameter:
        %         case {  'Fourier set','Fourier set (Hanning)',...
        %         'Gamma functions','Finite Impulse Response'}
        % Order, for FIR, is number of regressors to cover length period
        % order 6, length 30 -> 6 x 5 sec periods
        
        % bf = spm_get_bf(struct('name', 'hrf (with time and dispersion derivatives)', 'length', 30, 'dt', 1));
        
        bf = spm_get_bf(struct('name', varargin{2}, 'length', 30, 'dt', 1/res, 'order', 6));
        hrf = bf.bf;
        
        if doampscale && ~strcmp(varargin{2}, 'Finite Impulse Response')
            hrf = hrf ./ max(hrf);
        end
    
        % Eliminate varargin, so it won't interfere with getPredictors
        % later
        varargin{2} = '';
    else
        % custom HRF; do not norm by max
        % assume sampled in sec, and sample at high res first
        % scale to sum to 1 like SPM does
        hrf = varargin{2};
        lenh = length(hrf);
        hrf = interp1((1:lenh)', hrf, (1:(1/res):lenh), 'linear', 'extrap');
        hrf = hrf ./ sum(hrf);
    end
else
    
    hrf = spm_hrf(1/res, hrfparams);   % canonical hrf
    
    if doampscale
        hrf = hrf ./ max(hrf);
    end
end

% ----------------------------------------------
% Normalize basis set, if requested
% ----------------------------------------------

if donorm
    % Mean-center
    hrf = scale(hrf, 1);
    kk = size(hrf, 2);
    
    % Orthogonalize "later" basis funtions wrt "earlier" (assumed to be more important) ones
    if size(hrf, 2) > 1, hrf = spm_orth(hrf); end
    
    % Normalize to L2 norm
    for ii = 1:kk
        hrf(:, ii) = hrf(:, ii) ./ norm(hrf(:, ii));
    end
end

% ----------------------------------------------
% Set up custom neural response function
% ----------------------------------------------

if ~isempty(customneural)
    % assume sampled in sec, and sample at high res first
    % scale to sum to 1 like SPM does
    
    nsamp = length(customneural);
    customneural = interp1((1:nsamp)', customneural, (1:(1/res):nsamp), 'linear', 'extrap');
    
    if isrow(customneural), customneural = customneural'; end
    
    nsamp = length(customneural);
    
end

% ----------------------------------------------
% BUILD DESIGN
% see also getDesign5.m
% ----------------------------------------------

for i = 1:n_conditions

    cf2(:,i) = cf;                      % add empty column for this condition
    
    if isempty(customneural)
        % Event-related response (or epoch)
        
        cf2(round(ons{i}(:,1)*res) + 1, i) = 1;   % specify indicators; first TR is time 0, element 1
        
    else
        
        % Custom neural response
        for j = 1:n_events(i)
            
            first_sample = round(ons{i}(j, 1)*res) + 1;
            end_in_samples = min(len, first_sample + nsamp - 1);
            
            prevneural = cf2(first_sample:end_in_samples, i);
            cf2(first_sample:end_in_samples, i) = prevneural + customneural(1:length(prevneural)); % add to condition function
        
        end
    end
        
    if ons_includes_durations && ~isempty(customneural)
        error('You cannot specify epoch durations and a custom neural function.');
    end
    
    if ons_includes_durations
        % Build epochs
        % ----------------------------------------------------------------
        
        for j = 1:n_events(i)
            
            dur_in_s = ons{i}(j, 2);
            
            first_sample = round(ons{i}(j, 1) * res) + 1;
            end_in_samples = min(len, round(first_sample + dur_in_s * res));
            
            n_samples = dur_in_s .* res;
            
            % unmodulated: implicit assumption is that epoch is a series of
            % events of res (e.g., TR / 16. 
            % this gives very unrealistic relative amplitudes for brief vs.
            % long epochs. Within a trial type where this is constant, it
            % doesn't matter because scaling doesn't matter, but for
            % comparing across different event types with
            % onsets/epochs/blocks of different durations, and for
            % event-related designs with very rapid event densities (< 1
            % sec) then getting the relative scaling right is essential.
            % This influences simulations of efficiency, particularly at
            % very low ISIs (rapid events). One must implicitly assume an
            % event density for a block...what is the neural "frequency" of
            % events...1 per 1 msec? 1 per 1 sec?  
            %cf2(first_sample:end_in_samples, i) = 1; % add to condition function
            
            % modulated: equal area ('neural' function sums to 1)
            % ----------------------------------------------------
            % This is also very unrealistic, as it spreads one neural
            % "event" out over many seconds for long blocks.
            % The effective event density may vary from expt to expt
            % Nonlinear effects influencing the relative response amplitude
            % in Wager et al. 2005 are one reasonable starting point.
            % For designs with varying epoch durations, the effective event
            % density could be fit empirically so that the same linear
            % model captures both brief and long events.  This also applies
            % to parametric modulators, in which the change in relative
            % amplitude and shape with increasing neural frequency or
            % duration may often have nonlinear effects.
            %cf2(first_sample:end_in_samples, i) = 1 ./ n_samples; % add to condition function

            % modulated custom density: 
            % This assumes that the max neural sampling rate is about 1 Hz,
            % which produces responses that do not exceed about 5x the
            % unit, single-event response.  This is a reasonable
            % assumption, and affects only the absolute scaling across
            % ISIs/block lengths. 
            % ----------------------------------------------------
            cf2(first_sample:end_in_samples, i) = 1 ./ res; % add to condition function

        end
        
    end
    
    cf2 = cf2(1:len,:);
end

delta_hires = cf2;

% convolve and downsample. this now does parametric modulators as well.
% ------------------------------------------------------------------------
if ~isempty(len_original)
    % Downsample using 'dslen' input
    X = getPredictors(delta_hires, hrf, 'dslen', len_original/TR, 'force_delta', varargin{:}); % added len_original by Wani
else
    % Downsample using 'dsrate'
    X = getPredictors(delta_hires, hrf, 'dsrate', res*TR, 'force_delta', varargin{:}); % added len_original by Wani
end

if dononlinsaturation
    X = hrf_saturation(X);
end

X(:,end+1) = 1;


% for HRF shape estimation
% ------------------------------------------------------------------------

len2 = ceil(max(vertcat(ons{:})) ./ TR);
cf = zeros(len2(1), 1);
cf2 = [];

for i = 1:length(ons)
    tmp = 1 + round(ons{i}(:,1)./TR);
    cf2{i} = cf;
    cf2{i}(tmp) = 1;  % time 0 is first element
    
    repeats = tmp(find(diff(tmp) == 0));
    cf2{i}(repeats) = cf2{i}(repeats) + 1;
end

delta = cf2;

if ~isempty(varargin)
    for i = 1:length(ons)
        
        try
            % 2019 internal matlab function
            delta{i} = padarray(delta{i}, ceil(len./res./TR - length(delta{i})));
        catch
            % This is in CanlabCore/Misc_utilities
            delta{i} = pad(delta{i}, ceil(len./res./TR - length(delta{i})));
        end
        
    end
end

end % Main function


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Sub-functions

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


function [ons, n_conditions, n_events, ons_includes_durations] = check_onsets(ons, docheckorientation)
%
% [ons, n_conditions, n_events, ons_includes_durations] = check_onsets(ons, docheckorientation)
% n_conditions = number of event types, length(ons)
% n_events = vector of number of events in each condition

    % Enforce cell array
    
    if ~iscell(ons), ons = {ons}; end


    % sizes of onsets in each cell: First col is # events, 2nd is durations
    % if they exist.
    sz = cellfun(@size, ons, 'UniformOutput', false)';
    sz = cat(1, sz{:});
    
    % Check for empty cells / no events
    if size(sz, 1) < length(ons) || any(sz(:, 1) == 0)
        error('Some onsets are empty. Check inputs.');
    end
    
    % Check orientation
    % --------------------------------------------------
    if docheckorientation && all(sz(:, 1) == 1)
        
        % we have likely entered row vectors. transpose all.
        for i = 1:length(ons), ons{i} = ons{i}'; end
        
        sz = cellfun(@size, ons, 'UniformOutput', false)';
        sz = cat(1, sz{:});
    
    end

    n_conditions = size(sz, 1);
    n_events = sz(:, 1);
    
    % Check for illegal onsets
    % --------------------------------------------------
    for i = 1:n_conditions
        
    if any(round(ons{i}(:, 1)) < 0)
        error('Illegal onsets (negative values) -- check onsets!');
    end
    
    end
    
    ons_includes_durations = all(sz(:, 2) > 1);

    
end
