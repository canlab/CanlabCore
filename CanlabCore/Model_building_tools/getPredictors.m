function [model,delta] = getPredictors(stimList, HRF, varargin)
% Build predictors and delta functions, given a condition function or delta
% function and either a convolution matrix or vector.
%
% :Usage:
% ::
%
%     [model,delta] = getPredictors(stimList, HRF, varargin)
%
% IMPORTANT: YOU MUST ADD THE INTERCEPT YOURSELF!
%
% :Inputs:
%
%   **stimList:**
%        condition function OR delta function (1/0 indicator)
%
%   **HRF:**
%        1) hemodynamic response function
%        2) Basis set (columns)
%        3) or convolution matrix (columns
%           are HRF), defined as:
%           HRF = tril(toeplitz(hrf));
%
%           multiple column vectors for HRF are treated as basis functions!
%
%   **varargin** for downsampleing:
%        'dsrate': takes every nth element of the design matrix
%
%        'dslen': the target number (length) you want to downsample to
%
% :Other Optional Inputs:
%
%   **'force_delta':**
%        getPredictors tries to determine if the input stimList
%        is a condition function with integers or a delta function with
%        indicators, but this can fail in some cases.  Use this to force
%        it to treat as a delta function.
%
%   1. a col. vector of stimulus conditions OR a delta function matrix
%
%   2. an HRF vector sampled at the frequency of the stimulus vector, OR
%        a convolution matrix H (empty for default)
%
%   3. Optional: downsampling factor for final design (i.e., TR)
%
%   4. Optional: parametric modulator keyword and modulator values
%
%      'parametric_singleregressor' : Parametrically modulate onsets by
%               modulator values, using single regressor with modulated amplitude
%               Enter a cell array with modulator values for each event
%               type, with a column vector (empty cell for no modulation)
%
%       'parametric_standard' : Parametrically modulate onsets by
%               modulator values, using two regressors per event type - One
%               to model the average response, and one for the
%               mean-centered modulator values
%
% :Outputs:
%
%   1. a n x 2 matrix of regressors (cols) for each condition
%   2. a n x k delta matrix with onsets
%
% :Example: TR = 2, 16 samples per second in hi-res delta dhr
% ::
%
%    X = getPredictors(dhr, hrf, 'dsrate', res*TR); 
%    X = getPredictors(dhr, hrf, 'dslen', len/(res*TR)); 
%
% stimList can be condition function e.g., [1 3 2 4 3 2 1]' or
% delta matrix (n x k), n samples and k conditions, e.g., [1 0 0 0 1 0 1]'
%
% :Resampling: the default N in matlab resample has built-in antialiasing,
% but may not be good for fmri designs!  The appropriate downsampling
% is expected to be res*TR (res is units of samples/s), but we use 0
% because the model will depend on the analysis method used, and this is
% the most veridical approach.  With N = 0, every ith sample is used, where
% i is the downsampling factor you input.  Popular choices are 16*TR (for
% onsets2delta.m), using the SPM default res of 16.
% Delta is NOT resampled.
%
% :Example: TR = 2, 16 samples per second in hi-res delta dhr
% ::
%
%    [tmp,d] = downsample_delta(dhr,16*2); X=getPredictors(d,hrf);
%
% ..
%    Notes:
%    Tor Wager, last modified 2/22/04 to center predictors
%    modified to optionally take H convolution matrix as input
%    in place of HRF vector.  See convmtx.m (Buracas,mseq toolbox)
%
%    Modified 10/2011 by Tor to add parametric modulators
%
%    3/15/12: Tor updated to consider any stimlist with > 30 conditions to be
%    a delta matrix, so that one can input continuous neural response
%    functions and treat them as delta matrices for direct convolution.
% 
%    9/3/12: Wani updated this function for two reasons. 1) to fix an error 
%    when one condition has two or more parametric modulators. 2) to fix an
%    error due to RT with decimal points. For 2), make two options for downsampling
%   'dsrate' and 'dslen'.
% ..

model = [];
issquare = all(size(HRF) == size(HRF,1));   % if square mtx, assume convolution mtx

% make sure HRF is column vector, if single vector
if size(HRF,1) == 1, HRF = HRF'; end

% parametric modulators
doing_parametric_standard = any(strcmp(varargin, 'parametric_standard'));
doing_parametric_singleregressor = any(strcmp(varargin, 'parametric_singleregressor'));

if doing_parametric_standard || doing_parametric_singleregressor
    wh = find(strcmp(varargin, 'parametric_standard') | strcmp(varargin, 'parametric_singleregressor'));
    pm_vals = varargin{wh(1) + 1};
    
    if ~iscell(pm_vals)
        error('pm_vals must be cell array.')
    end
end

use_downsample_rate = any(strcmp(varargin, 'dsrate'));
use_downsample_length = any(strcmp(varargin, 'dslen'));
if use_downsample_rate, wh = find(strcmp(varargin, 'dsrate')); rate_downsample = varargin{wh(1) +1}; end
if use_downsample_length, wh = find(strcmp(varargin, 'dslen')); length_downsample = varargin{wh(1) +1}; end


if isdeltamtx(stimList, varargin{:})  % delta matrix
    % -------------------------------------------------------------------------------------------------
    %    * If delta function
    % -------------------------------------------------------------------------------------------------
    
    delta = stimList;
    
    if ~issquare
        for i = 1:size(delta,2)
            
            for j = 1:size(HRF,2)
                
                if doing_parametric_standard || doing_parametric_singleregressor
                    if length(pm_vals) < i, error('PM values cell array is too short!'); end
                    
                    model = [model pmconv(double(delta(:,i)), HRF(:,j), pm_vals{i}, doing_parametric_standard, doing_parametric_singleregressor)];
                else
                    model(:,end+1) = conv(double(delta(:,i)), HRF(:,j)); % Changed for Matlab 7.9+ compatibility - Thanks, Liane
                end
            end
        end
    end
    
else
    % -------------------------------------------------------------------------------------------------
    %    * If condition function
    % -------------------------------------------------------------------------------------------------
    
    for i = 1:max(stimList(:,1)) % condition function
        
        %delta(:,i) = (stimList == i);
        delta(:,i) = cast((stimList == i), 'single'); % Changed for Matlab 7.9+ compatibility -- thanks, Liane Montana and Bruce McCandliss
        
        if doing_parametric_standard || doing_parametric_singleregressor
            if length(pm_vals) < i, error('PM values cell array is too short!'); end
            
            model = [model pmconv(double(delta(:,i)), HRF(:,j), pm_vals{i}, doing_parametric_standard, doing_parametric_singleregressor)];
        end
        
        if ~issquare
            for j = 1:size(HRF,2)
                model(:,end+1) = conv(double(delta(:,i)), HRF(:,j)); % Changed for Matlab 7.9+ compatibility - Thanks, Liane
            end
        end
        
    end
end

% -------------------------------------------------------------------------------------------------
%    * If conv. matrix
% -------------------------------------------------------------------------------------------------

if issquare % convolution matrix
    model = HRF * delta;
end

model = model(1:size(stimList,1),:);      	% eliminate extra values

% downsample, if necessary
if use_downsample_rate || use_downsample_length % Wani modified this
% if ~isempty(varargin)    
    
    % dsrate = varargin{1}; % Wani
    [n, k] = size(model); 
    % nt = n ./ dsrate; % Wani
    if use_downsample_rate, nt = n ./ rate_downsample; end % Wani modified this line
    if use_downsample_length, nt = length_downsample; rate_downsample = n/length_downsample; end % Wani added this line
    t = (1:n)';             % minor change to save time, tor: 10/12/10
    
    if nt-round(nt)>.000001, error('Length of stimList is not evenly divisible by downsampling factor.'); end
    
    modeli = zeros(round(nt) , k);
    
    for i = 1:size(model, 2)
        
        xi = 1:rate_downsample:n;       % downsample rate
        
        modeli(:, i) = interp1(t, model(:, i), xi, 'linear', 'extrap'); % tor: 10/12/10 : minor change: add extrap
    end
    
    model = modeli;
    %model = model(1:varargin{1}:end,:); % equivalent to resample(X,1,varargin{1},0)
end

% do not do this if you want to do nonlinear saturation!
% not necessary, i think.
%model = model - repmat(mean(model),size(model,1),1);    % center predictors

end




function isdelta = isdeltamtx(stimList, varargin)

isdelta = 0;

if ~isempty(varargin) && any(strcmp(varargin, 'force_delta'))
    isdelta = 1;
    return
end


if islogical(stimList) || all( sum(stimList == 0 | stimList == 1) == size(stimList, 1) )  % all zeros or ones
    isdelta = 1;
    
elseif min(size(stimList)) > 1 % multiple columns
    isdelta = 1;
    
elseif length(unique(stimList)) > 30 % custom neural response function?
    
end

end % subfunction



function [model, delta] = pmconv(delta, HRF, pm_vals, doing_parametric_standard, doing_parametric_singleregressor)

if doing_parametric_standard && doing_parametric_singleregressor
    error('Cannot do both standard and single-regressor parametric modulator.')
end

% HRF can be basis set, but SINGLE basis function is entered
% this is for one regressor/event type, with multiple basis images

% must identify FIRST events, i.e., preceded by 0
% for epochs... and must get back epoch duration, too.
%wh = find(delta);
wh = find(delta & [1; diff(delta) == 1]);

% wh_end values, for durs
wh_end = find(~delta & [0; diff(delta) == -1]);
if delta(end), wh_end(end+1) = length(delta); end % fix if last event goes to end

% transpose if horiz
if diff(size(pm_vals)) > 0, pm_vals = pm_vals'; end

if ~isempty(pm_vals) && length(wh) ~= length(pm_vals)
    disp('Number of events in parametric modulator does not match number of events in delta function.')
    fprintf('%3.0f events in delta function, %3.0f events in param. modulator\n', sum(delta), length(pm_vals));
    
    disp('Bad input?  Stopping in debugger so you can check. Type dbcont or dbquit to continue.')
    keyboard
end

% get modulated delta function

if ~isempty(pm_vals) && doing_parametric_standard
    % create a second delta function that is modulated
    pm_vals = scale(pm_vals, 1);
    delta = repmat(delta, 1, size(pm_vals,2)+1); % Wani modified this line to fix an error when pm_vals has two or more columns. 
    % delta(:, end + 1) = delta; 
    
elseif ~isempty(pm_vals) && doing_parametric_singleregressor
    % do nothing; we will just use pm vals
    
elseif isempty(pm_vals)
else
    error('This should not happen');
end

% This is where they get added.
% tor adjusted aug 2012 to account for epoch case.
if ~isempty(pm_vals)
    for i = 1:length(wh)
        delta(wh(i):wh_end(i), end) = pm_vals(i);
    end
end

% Convolve all preds for this trial type
lendelta = size(delta, 1);

for i = 1:size(delta, 2)
    
    %for j = 1:size(HRF, 2)
    tmp = conv(double(delta(:,i)), HRF); % Changed for Matlab 7.9+ compatibility - Thanks, Liane
    
    X{i} = tmp(1:lendelta);
    %end
    
end

model = cat(2, X{:});


end


