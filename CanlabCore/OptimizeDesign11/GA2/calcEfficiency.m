function [eff,eff_vector, contrasts] = calcEfficiency(contrastweights,contrasts,xtxitx,svi,varargin)
% Calculate the weighted average of the design efficiency values for a set of contrasts.
% :Usage:
% ::
%
%     [eff,eff_vector, contrasts] = calcEfficiency(contrastweights,contrasts,xtxitx,svi,[D-opt flag])
%
% - Used in GA toolbox
% - Inputs configured for speed in repeated tests and backwards compatibility with GA toolbox
% - Can calculate average efficiency for individual regressors or contrasts
% - Can automatically create a set of orthogonal contrasts spanning the space of differences among events (conditions)
%
% - For contrast matrix C (num_cons x num_regressors), efficiency for each contrast is defined as: 
%   1 ./ ( ( diag(C * pinv(X) * svi * pinv(X)' * C') ./ diag(contrasts * contrasts') ) .^ .5 )
%
% :Inputs:
%
%   **Fixed set of inputs in order:**
%
% contrastweights   a row vector of length equal to number of rows of contrasts indicating relative weights
%                   If a vector of ones, this will calculate the average
%                   contrasts (the function divides by n contrasts)
%                   If empty, this function will create a vector of 1/n_contrasts
% contrasts         a matrix of contrasts to test, one row per contrast, one col per condition in X
% xtxitx            pinv(X), where X is the design matrix
% svi               S * Vi * S', where S is the filtering matrix and Vi is the intrinsic autocorrelation matrix
%                   Can be empty, in which case a default matrix is used based on AR(1)
%
% eff               overall efficiency, equal to eff_vector * contrastweights for A-optimality
%                   or det[Z] for D-optimality, where Z is a matrix you can see if you look at the code
% eff_vector        vector of efficency scores for each contrast
% [D-opt flag]      optional argument that specifies computation of D-optimality,
%                   which is the det of the inversion of the Fisher Info Matrix instead of the trace
%                   this is better for model selection, whereas A-opt is better for testing a particular
%                   contrast
%                   % IF USING THIS FLAG, xtxitx SHOULD BE THE FISCHER INFO MTX X'X !!!!
%                   THIS DOPT THING IS NOT WORKING PROPERLY, SO I WOULDN'T USE THE D-OPT FLAG HERE AT THIS TIME.
%                   see ga2.m for a simple way to compute d-opt.
%
%   **Alternative - single input X or 2 inputs:**
%
%   X is a design matrix with a single intercept.  The function:
%   - Enforces that the intercept is at the end
%
%   [eff,eff_vector, contrasts] = calcEfficiency(X);
%   - Calculates average efficiency for individual regressors
%   - The efficiency of the intercept is the last element of eff_vector
%   - Efficiency of the intercept is not used in calculation of eff
%
%   [eff,eff_vector, contrasts] = calcEfficiency(X, 'contrasts');
%   - Calculates average efficiency for orthogonal contrasts capturing
%   differences among regressors
% 
%   [eff,eff_vector, contrasts] = calcEfficiency(X, C);
%   - For contrast matrix C [n_contrasts x regressors including intercept]
%   - Calculates average efficiency for a priori contrasts
%
% by Tor Wager. Circa 2005, updated 2021.
%
% Examples:
% ------------------------------------------------------------------------
% Create a sparse event-related design with an event every 10 sec, TR =
% 1.3, events are 2-sec epochs:
% ons = mat2cell([1:10:200]', ones(20, 1))';
% [X,d,out,handles] = plotDesign(ons,[], 1.3, 'durs', 2);
%
% % Get efficiency (average efficiency for events)
% [eff,eff_vector, contrasts] = calcEfficiency(X);
%
% % Get efficiency (create orthogonal set of contrasts)
% [eff,eff_vector, contrasts] = calcEfficiency(X, 'contrasts');
%
% ------------------------------------------------------------------------
% % Create random event-related design with an event every 3 sec
% % 4 conditions with 20% frequency each, 2 sec epochs, and plot it with 1.3 sec TR:
% ons = create_random_onsets(100, 3, [.2 .2 .2 .2], 2);
% [X,d,out,handles] = plotDesign(ons,[], 1.3);
%
% % Plot VIFs for this design:
% create_figure('vifs'); getvif(X, 0, 'plot');
%
% % Get efficiency
% nconditions = length(ons);
%
% contrasts = create_orthogonal_contrast_set(nconditions);
% contrasts(:, end+1) = 0; % for intercept
%
% e = calcEfficiency(ones(1, size(contrasts, 1)), contrasts, pinv(X), []);
%
% ------------------------------------------------------------------------
% % Create a dense random event-related design, 2 sec epochs, 1.3 sec TR
% ons = create_random_onsets(100, 2, [.5 .5], 2);
% [X,d,out,handles] = plotDesign(ons,[], 1.3);
% 
% % Average efficiency for events (vs intercept) 
% [eff,eff_vector, contrasts] = calcEfficiency(X); eff
% 
% % Efficiency for auto-generated contrast set for diffs among conditions (here, just [1 -1 0])
% [eff,eff_vector, contrasts] = calcEfficiency(X, 'contrasts'); eff
% 
% % Average efficiency for diff and sum of two events
% [eff,eff_vector1, contrasts] = calcEfficiency(X, [1 -1 0; 1 1 0]); eff
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2005 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..

% Programmers' notes:
% Documentation update: 9/18/09
% Code cleanup, one minor change to estimate when no contrasts are entered
% (use n - 1) because we assume an intercept.
%
% Note: efficiency depends on scaling of contrastweights, contrasts, and
% hrf.  Comparing different models (e.g., with different basis sets) may be very tricky.
% in this implementation, contrastweights should sum to
% length(contrastweights) and be all positive for comparing designs.
%
% Power is proportional to sqrt(efficiency), holding constant the:
%  1) "optimal" effect size, magnitude over sigma
%  2) group sample size, and assuming between-subjects error = 0
% t = (M/sigma)/(1/sqrt(var)) =  (M/sigma) * sqrt(eff) <- if eff is in std units then it's just eff
% assume M/sigma = 1, t = sqrt(eff) <- if eff is in std units then it's just eff
% follows non-central t distribution with mean M/sigma
% p = 1 - tcdf(Z, df)
% p = 1 - tcdf(sqrt(eff), size(X, 1) - size(X, 2))
% -log(p) is -log of expected p-value given design and std. effect size
%
% Dec 2021: Tor Wager:
% Added input options to simplify
% Changed the scaling of efficiency in two ways:
% (1) Divides by the inner product of contrast weights(C * C') so
% efficiency is invariant wrt contrast weights
% (2) Expressed in units of standard dev (sqrt(var)) rather than variance

Dflag = 0;

% ------------------------------------------------------------------------
% Special case with X and possibly contrasts only
% ------------------------------------------------------------------------
if nargin == 1 || nargin == 2
    
    % Special case where design matrix only is entered
    % Argument structure is still fixed for backwards compatibility
    
    X = contrastweights;        % the first input
    X = intercept(X, 'end');    % enforce intercept as last regressor

    % Uses default autocorrelation and smoothing, assuming TR = 1 and 180 sec HP filter
    % These choices will affect the efficiency of low-frequency designs,
    % and even if the TR varies from this it will often be a reasonable
    % approximation.
    
    [S, ~, svi] = getSmoothing(180, [], 1, size(X, 1), 'auto');
    X = S*X;
    
    % Alternative to the two lines above:
    % svi = [];
     
    xtxitx = pinv(X);
    contrastweights = [];
    
    if nargin == 1
        % one input; empty contrasts
        contrasts = [];
        
    else
        % two inputs: we have entered something for contrasts...parse this
        % Contrast weights, string 'contrasts', or something else (ignore
        % the latter case)
        
        if ~ischar(contrasts)
            % assume these are contrast weights
            % do nothing now
            
        else
            % Enter 'contrasts' to estimate average over a set of
            % orthogonal contrasts, or no 2nd argument to estimate average
            % for each condition.
            
            % if string, create contrasts
            if strcmp(contrasts, 'contrasts')
                
                nconditions = size(X, 2) - 1;
                
                contrasts = create_orthogonal_contrast_set(nconditions);
                contrasts(:, end+1) = 0; % for intercept
                
                %                 contrastweights = ones(1, size(contrasts, 1)); % average contrast
                
            else
                % We have 2nd argument, but empty, so we use average condition
                contrasts = [];
                %                 contrastweights = [];
                
            end % String contrast input or nothing
            
        end % contrast handling
        
    end % 1 or 2 inputs
    
end % Special input case

% ------------------------------------------------------------------------
% Proceed with Standard case with fixed inputs
% ------------------------------------------------------------------------

if isempty(svi), svi = eye(size(xtxitx,2)); end

if length(varargin) > 0, Dflag = varargin{1}; end

% if empty use a vector of ones for the number of contrasts
if isempty(contrastweights), contrastweights = ones(1, size(contrasts, 1));  end

if ~isempty(contrasts)
    % compute the variance of the contrasts
    % ---------------------------------------------------------------------
    
    try
        % Variances of each contrast: (note: reciprocal taken later)
        % ---------------------------------------------------------------------
        scaling_factor = diag(contrasts * contrasts');
        var_vector = diag(contrasts * xtxitx * svi * xtxitx' * contrasts') ./ scaling_factor;
  
    catch
        disp('Error calculating efficiency: design matrix, filtering matrix, and contrasts matrix sizes do not match!!!')
        disp('Trying to do this calculation: eff_vector = diag(contrasts*xtxitx*svi*xtxitx''*contrasts'')')
        
        contrastsize = size(contrasts);
        pxsize = size(xtxitx);
        svisize = size(svi);
        
        disp(' ');
        disp('I am trying to multiply these matrices:')
        fprintf('C (%01d x %01d) * PX (%01d x %01d) * SVi (%01d x %01d)\n ', contrastsize, pxsize, svisize)
        
        disp('The error is in this multiplication.')
        disp('If any adjacent numbers do not match up in the line above, something has been entered incorrectly')
        disp('The most common cause is a contrast matrix that has either too many or too few columns.')
        error('Exiting.');
    end
    
    try
        if Dflag
            eff = 1 ./ inv(det(xtxitx)) .^ (1./size(xtxitx, 2));   % HERE xtxitx SHOULD BE THE FISCHER INFO MTX X'X
        else
            % compute weighted average over contrasts depending on
            % contrastweights. Take reciprocal of variance, so instead of wavg/n
            % use n/wavg.  If contrastweights are all ones, will
            % average.
            % ---------------------------------------------------------------------
            
            eff = length(var_vector) ./ (contrastweights * var_vector);
        end
        
        
    catch
        disp('Error calculating efficiency: contrast weights and number of contrasts do not match!!!')
        disp(' I am trying to do this calculation:')
        disp('eff = length(eff_vector)./(contrastweights * eff_vector);')
        
        cwsize = size(contrastweights);
        effsize = size(var_vector);
        
        disp(' ');
        disp('I am trying to multiply these:')
        fprintf('contrastweights (cwsize) * eff_vector (%01d x %01d) \n ', cwsize, effsize)
        disp('The most common causes are too many or too few entries in contrastweights.');
        error('Exiting.');
    end
    
else
    
    % contrasts are not entered.  Assume we want the average efficiency of
    % each individual model parameter.
    %
    % ---------------------------------------------------------------------
    
    % If contrasts are empty, contrastweights may still be empty:
    if isempty(contrastweights)
        contrastweights = [ones(1, size(xtxitx, 1) - 1) 0]; % one for each beta, including intercept
    end
    
    try
        % Variances of each regression parameter:
        % ---------------------------------------------------------------------
        var_vector = diag(xtxitx * svi * xtxitx');
        
        if Dflag
            
            eff = 1 ./ inv(det(xtxitx)) .^ (1./size(xtxitx, 2));   % HERE xtxitx SHOULD BE THE FISCHER INFO MTX X'X
            
        else
            % compute the average variance of the chosen regressors
            % Take reciprocal so higher values = more efficient
            % Use n - 1 because we are assuming an intercept, and have one
            % irrelevant parameter.
            % ---------------------------------------------------------------------
            eff = (length(var_vector) - 1) ./ (contrastweights * var_vector);
            
        end
        
        
    catch
        disp(['no contrasts version: wrong sizes!!!'])
        %contrastweights = [contrastweights 0]
        whos contrastweights
        whos xtxitx
        whos svi
        diag(xtxitx*svi*xtxitx')
        error('Exiting...')
    end
    
end % contrasts yes/no

% Scale to standard deviation rather than variance
% This is directly proportional to the t-statistic then
eff = eff .^ .5;

% Take reiprocal of individual variances
% because we want to minimize the variance, but we maximize
% fitnessMatrix value.
eff_vector = 1 ./ (var_vector .^ .5);


end % function