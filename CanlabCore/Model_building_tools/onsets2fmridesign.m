% [X, delta, delta_hires, hrf] = onsets2fmridesign(onsets, TR, [len], [custom hrf or basis set name],[other optional args, i.e., norm])
% Tor Wager, 2 / 24 / 04
% Modified/updated 2/08 by Tor to add capacity for durations
% Modified 10/1/11 by Tor to fix custom HRF sampling
% Modified 3/14/12 by Tor to check variable duration and add noampscale option
%
% Builds design matrix X and delta function, given cell array of onset times in s
% One cell per condition per session, e.g., ons{1} = [24 27 29 44]';
%
% Summary:
% - handles multiple conditions
% - handles custom HRFs and multiple basis functions
% - handles two kinds of parametric modulators
% - handles variable-duration onsets
% - not yet: variable-duration parametric
%       modulators
%
%
% Inputs:
%   onsets:
%   - 1st column is onsets for events,
%   - 2nd column is optional durations for each event
%   - Enter single condition or cell vector with cells for each condition (each event type).
%   TR: RT in seconds
%   len: optional length in s for model, or [] to skip if using additional. 
%        "len" is usually the number of images multiplied by TR.
%
%   args
%   hrf name: 1) a string used by spm_get_bf.m
%             2) a custom HRF, sampled in seconds
%             3) or [] to skip
%   'norm' : mean-center, orthogonalize, and L2-norm basis set
%   'parametric_singleregressor' : Parametrically modulate onsets by
%               modulator values, using single regressor with modulated amplitude
%               Enter a cell array with modulator values for each event
%               type, with a column vector (empty cell for no modulation)
%   'parametric_standard' : Parametrically modulate onsets by
%               modulator values, using two regressors per event type - One
%               to model the average response, and one for the
%               mean-centered modulator values
%               of modulator values in each cell
%   'noampscale' Do not scale HRF to max amplitude = 1; SPM default is not to scale.
%   'noundershoot' Do not model undershoot - ONLY when using the canonical HRF
%   'customneural' Followed by a vector of custom neural values (instead of
%                   standard event/epochs), sampled in sec.
%
% Limitation: can only handle two events max within the same TR
%
% Outputs:
%   X: model, sampled at TR
%   delta: indicator mtx, sampled at TR
%   delta_hires: indicator sampled at 16 Hz
%   hrf: hemodynamic response function, sampled at 16 Hz
%
% Test examples:
% X = onsets2fmridesign(ons, TR, size(imgs, 1) .* TR, 'hrf (with time derivative)');
%
% X = onsets2fmridesign({[0 30 60]' [15 45 90']}, 1.5, 180, spm_hrf(1), 'parametric_standard', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
% X = onsets2fmridesign({[0 30 60]' [15 45 90']}, 1, 180, spm_hrf(1), 'parametric_standard', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
% X = onsets2fmridesign({[0 30 60]' [15 45 90']}, 2, 180, spm_hrf(1), 'parametric_standard', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
% X = onsets2fmridesign({[0 30 60]' [15 45 90']}, 2, 180, spm_hrf(1), 'parametric_singleregressor', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
% X = onsets2fmridesign({[0 30 60]' [15 45 90']}, 1, 180, spm_hrf(1), 'parametric_singleregressor', {[2 .5 1]' [1 2 2.5]'}); figure; plot(X)
%
%   Here, spm_hrf(1) is for the canonical HRF in spm.
%
% X = onsets2fmridesign(ons, 2, size(dat.dat, 2) .* 2, [], [], 'noampscale');
%
% X2 = onsets2fmridesign(onsp, 2, length(dat.removed_images) .* 2, [], [], 'customneural', report_predictor);
%
% Programmers' notes:
% 3/16/12: Tor modified custom HRF to divide by sum, as SPM basis functions
% are scaled by default.
% 9/3/12: Wani modified len and getPredictors.m to fix a bug when you're
% using TR with some decimal points (e.g., TR = 1.3)

function [X, delta, delta_hires, hrf] = onsets2fmridesign(ons, TR, varargin)
% ----------------------------------------------
% Defaults
% ----------------------------------------------

res = 16;   % resolution, in samples per second

% Other optional inputs after fixed ones
% - - - - - - - - - - - - - - - - - - - -
donorm = 0;
doampscale = 1;
hrfparams = [6 16 1 1 6 0 32];
customneural = []; % custom neural response; neither impulse nor epoch
dosingletrial = 0;
docheckorientation = 1;

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
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


if ~iscell(ons), ons = {ons}; end

if dosingletrial
    % single trial: break events into separate regressors/onset vectors
    for i = 1:length(ons)
        sonsets{i} = cell(1, length(ons{i})); % event-specific onsets - init
        
        for j = 1:length(sonsets{i})
            sonsets{i}{j} = ons{i}(j, 1);
            if size(ons{i}, 2) > 1   % if durations...
                sonsets{i}{j} = [sonsets{i}{j} ons{i}(j, 2)]; % append duration
            end
        end
        
    end % original onset conditions
    ons = cat(2, sonsets{:}); % flatten to simply model trials within event type.
end


% Optional: pre-specified length, or []
% - - - - - - - - - - - - - - - - - - - -

if ~isempty(varargin) && ~isempty(varargin{1}) % pre-specified length
    
    len_original = varargin{1}; % this line is added by Wani to make this work for TR = 1.3
    len = varargin{1}; % this line is modified by Wani from the next line
    % len = ceil(varargin{1});
    for i = 1:length(ons)
        if docheckorientation
            % make column vectors if needed
            if length(ons{i}) > size(ons{i}, 1)
                disp('Warning! onsets2fmridesign thinks you''ve entered row vectors for onsets and is transposing.  Enter col vectors!');
                ons{i} = ons{i}';
            end
        end
        
        % make sure nothing goes past end of session
        past_end = round(ons{i}(:,1)) > len;
        if any(past_end)
            fprintf('Warning! onsets2fmridesign found %d onsets after the end of the session and is omitting them.\n', sum(past_end));
        end
        ons{i}(past_end,:) = [];
    end
    len = ceil(len .* res);
else
    len = round(max(cat(1, ons{:})) * res);
    len = ceil(ceil(len/(res*TR)) * res * TR); % modified by Wani to make this work for TR = 1.3
    len = len(1);
end

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

for i = 1:length(ons)
    
    if any(round(ons{i}(:,1)) < 0)
        error('Illegal onsets (negative values) -- check onsets!');
    end
    
    cf2(:,i) = cf;                      % add empty column for this condition
    
    if isempty(customneural)
        % Event-related response (or epoch)
        cf2(round(ons{i}(:,1)*res) + 1, i) = 1;   % specify indicators; first TR is time 0, element 1
        
    else
        % Custom neural response
        for j = 1:size(ons{i}, 1)
            
            first_sample = round(ons{i}(j, 1)*res) + 1;
            end_in_samples = min(len, first_sample + nsamp - 1);
            
            prevneural = cf2(first_sample:end_in_samples, i);
            cf2(first_sample:end_in_samples, i) = prevneural + customneural(1:length(prevneural)); % add to condition function
        end
    end
    
    we_have_durations = size(ons{i}, 2) > 1;
    
    if we_have_durations && ~isempty(customneural)
        error('You cannot specify epoch durations and a custom neural function.');
    end
    
    if we_have_durations
        % Build epochs
        for j = 1:size(ons{i}, 1)
            dur_in_s = ons{i}(j, 2);
            first_sample = round(ons{i}(j, 1)*res) + 1;
            end_in_samples = min(len, round(first_sample + dur_in_s * res));
            
            cf2(first_sample:end_in_samples, i) = 1; % add to condition function
        end
        
    end
    
    cf2 = cf2(1:len,:);
end

delta_hires = cf2;

% convolve. this now does parametric modulators as well.
X = getPredictors(delta_hires, hrf, 'dslen', len_original/TR, varargin{:}); % added len_original by Wani
X(:,end+1) = 1;


% for HRF shape estimation
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
        delta{i} = pad(delta{i}, ceil(len./res./TR - length(delta{i})));
    end
end
end

