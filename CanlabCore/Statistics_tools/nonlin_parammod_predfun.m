function [yhat,yhi,loTime,highTime,x] = nonlin_parammod_predfun(y,ons,pm_vals,p,varargin)
% Predict data given onsets, parametric modulator values, and parameter
% estimates for effects of events and PMs on the magnitude and delay of an
% estimated neural 'boxcar.'
% 
% A standard HRF model is assumed, and the parameters are estimates of how
% trials affect the 'neural' stimulus function.
% 
% This function replaces nonlin_parammod_predfun.m, which is obsolete
%
% :Usage:
% ::
%
%     [yhat,yhi,loTime,highTime,x] = nonlin_parammod_predfun(y,ons,pm_vals,p,varargin)
%
% :Inputs:
%
%   **p:**
%        params
%          1) mag. scale (amplitude)
%          2) mag. X RT slope (linear effect of RT)
%          3) duration intercept (neural epoch duration)
%          4) duration x RT slope (linear effect of RT on duration)
%          5) overall fitted response intercept
%
%   **y:**
%        data. Can be dummy data, needs only be the correct size
%        This function does not actually fit the data, so y is only used
%        to get the correct vector size
%
%   **ons:**
%        onsets (in samples)
%
%   **pm_vals:**
%        modulator values, preferably centered
%
%   **p:**
%        parameter values; see below
%
%   **varargin:**
%        optional inputs, see below
%
% :Special cases of parameter sets:
%
% p = [1 0 1 0 0];  % "Impulse model" : fixed boxcar of 1 s
%
% p = [1 0 3 0 0];  % "Epoch model" : fixed boxcar of 3 s
%
% p = [1 1 1 0 0];  % "Parametric modulator" : Impulse height modulated by
% RT (with centered pm_vals)
%
% p = [1 0 0 1 0];  % "Variable epoch": convolve RT with duration 
% (with non-centered, raw pm_vals).  Parametric modulation of duration
%
% p = [1 1 1 1 0];  % combo of Parametric and Duration modulators of impulses
%
% p = [1 1.5 3 1.5 0];  % combo of Parametric and Duration modulators of
%                   a 3-s epoch model
%
% event signal magnitude = p1 * p2*RT           % was previously: RT^p2
% event signal  duration = p3 + pm_vals * p4
%
% :Optional Inputs:
%   - case 'random', dorandom = 1;
%   - case 'plot', doplot = 1;
%   - case 'hrf', hrf = varargin{i + 1}; % HIGH_RES HRF
%   - case 'tr', followed by TR
%
% :Examples of generating predicted BOLD timeseries:
% ::
%
%    y = rand(330,1);
%    ons = [1:20:320]'; pm_vals = rand(size(ons));
%    [yhat,yhi] = nonlin_parammod_predfun(y,ons,pm_vals,[1 1 1 0 0],'plot');  % linear RT modulation of height
%    [yhat,yhi] = nonlin_parammod_predfun(y,ons,pm_vals,[1 0 1 10 0],'plot');  % linear RT modulation of duration
%
%    [yhat,yhi] = nonlin_parammod_predfun(y,ons,pm_vals,[.5 1 1 0 0],'plot');  % other linear RT modulation of height
%    [yhat,yhi] = nonlin_parammod_predfun(y,ons,pm_vals,[.5 .7 1 5 0],'plot');  % saturated RT modulation of height; linear width
%
% Example of fitting observed timeseries and estimating parameters:
%
% * SEE ALSO nonlin_param_modulator
% ::
%
%    RTcenter = pm_vals - mean(pm_vals);
%    true_p = [1 1.5 3 1.5 0];          % true parameters
%    xvals = (1:length(y))';
%
% Fitting function: times is a dummy var to get this to work with nonlin_fit
% ::
%
%    fhan = @(p, times) nonlin_parammod_predfun(y,ons,RTcenter,p);
%
%    y = fhan(true_p);                  % generate simulated "true" signal
%
%    fitting_fun = @(y) nonlin_fit(y, xvals, 'link', fhan, 'start',[1 1 1 1 .5]);
%    [p, errval, fit] = fitting_fun(y)     % get parameters for a timeseries of interest
%
%    [p,errval,fit] = nonlin_fit(y,(1:length(y))','link',fhan,'start',[1 1 1 1 .5]);
%    hold on; plot(fit,'r');
%
% A second example using the genetic algorithm:
% ::
%
%    objfun = @(p) sum(abs(y - fhan(p))); % objective function: absolute error
%    start = [-10 -10 0 -10 0]; start(:,:,2) = [10 10 10 10 2000];
%    [best_params,fit,beff,in] = tor_ga(200,50,{start},objfun,'genconverge',5,'noverbose');
%    %***note: does not work now, needs debugging***
%
% Simulation: Test cov of param estimates
%
% Run 1000 times to see cov. of param estimates
% Right now: ****development: p(1,2) are highly + corr, p(3:4) are high - corr
% should probably choose one or the other for each.
% ::
%
%    for i = 1:1000
%       yi = y + randn(length(y),1) * 10; % new noise
%       p(i,:) = fitting_fun(yi);
%       if mod(i,10) == 0, fprintf(1,' %3.0f',i); end
%    end
%
% Note: Could create linear basis set of plausible forms
%
% ..
%    tor wager, jan 07
% ..

if ~isrow(pm_vals), pm_vals = pm_vals'; end

loResolution = 1;   % in seconds (this is the TR)
hiResolution =.1;  % in seconds
relSampRate = round(loResolution./hiResolution);
runDuration = 0;

% setup inputs
dorandom = 0;
doplot = 0;
if length(varargin) > 0
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'random', dorandom = 1;
                case 'plot', doplot = 1;

                case 'hrf', hrf = varargin{i + 1};
                    
                case 'tr'
                    loResolution = varargin{i + 1};
                    relSampRate = round(loResolution./hiResolution);
                    
            end
        end
    end
end

if dorandom
    setup_random;
else
    setup_vars;
end

% m, magnitude exponent
% p(2) = 2) mag. X RT exponent
% magnitude = p(1) .* pm_vals .^ p(2); % p(2) .^ pm_vals; %pm_vals .^ p(2);
% duration = relSampRate .* (p(3) + p(4) .* pm_vals); % in hi-res samples ; could be .* pm_vals, but then neg. durations are possible?

% want RT effects to be independent from overall trial effects
% with centered pm_vals, this should be accomplished by:

magnitude = p(1) + p(2) .* pm_vals;  % .^ p(2); 
duration = relSampRate .* (p(3) + p(4) .* pm_vals);

% but pm_vals must be positive to use exponents, e.g., pm_vals .^ p(2), otherwise fractional p's will
% give imaginary numbers.

% specify 'neural' boxcars based on mag and dur
xbox = get_boxcars(ons, x, magnitude, duration);

if length(hrf) > length(xbox), hrf = hrf(1:length(xbox)); end
    
% convolve
yhi = p(5) + fast_conv_fft(hrf, xbox);

if length(yhi) < runDuration
    yhi = pad(yhi, runDuration);
elseif length(yhi) > runDuration
    yhi = yhi(1:runDuration);
end

% downsample and add intercept
yhat = yhi(1:relSampRate:end);

n = length(y);
if length(yhat) > n, yhat = yhat(1:n); end

if doplot || nargout > 2
    loTime = (1:length(y))' .* loResolution - 1;
end

if doplot
    create_figure('Predicted fMRI Response');
    barht = max(yhi) ./ (4.*max(x));
   
    %tor_fig;
    plot(highTime,barht.*x,'k',highTime,yhi,'b',loTime,yhat,'r');
    bar(highTime,barht.*x,'FaceColor',[.4 .4 .4],'EdgeColor','none');
end

% END MAIN FUNCTION


    function setup_random
        runDuration=330;    % in seconds
        highTime=linspace(0, runDuration, runDuration/hiResolution)';

        hrf=gam(highTime,.2,.5,3,-.06,3,7,0);
        hrf=hrf/max(hrf);

        % get onsets
        numTrials = 20;
        x = zeros(3300,1);
        ons = randperm(length(x)); ons = ons(1:numTrials);
        x(ons) = 1;

        % % z = zeros(length(hrf),1); z(1) = hrf(1);
        % % HRF = toeplitz(z,hrf);

        % get pm_vals
        a=2.5;
        b=.13;
        minRT=.5;
        pm_vals=gamrnd(a,b,numTrials,1)+minRT;

    end


    function setup_vars

        runDuration = relSampRate .* length(y); % in elements
        x = zeros(runDuration,1);

        % ons entered in TRs, convert to hi-res
        % would multiply by 1./hiRate if entered in s
        ons = round(relSampRate .* ons);
        ons(ons == 0) = 1; % onsets at zero -> first high-res sample
        
        % setup HRF
        highTime = (0:runDuration-1)' .* hiResolution; %linspace(0, runDuration, runDuration/hiResolution)';

        %hrf=gam(highTime,.2,.5,3,-.06,3,7,0);
        if ~exist('hrf', 'var')
            hrf = spm_hrf(hiResolution);
            hrf = hrf/max(hrf);
        end
        
        hrf = pad(hrf,x);
        
        % Check
        if length(pm_vals) ~= length(ons)
            error('Length of pm_vals and length of ons must match.');
        end

    end


end  % END MAIN FUNCTION


function x = get_boxcars(ons,x,magnitude,duration)

duration = round(duration);

% x is high-res indicator

for i = 1:length(ons)
    % for each onset
    st = ons(i);
    en = ons(i) + duration(i) - 1;
    
    x(st : en) = magnitude(i);
end


end
