function [eff,eff_vector] = calcEfficiency(contrastweights,contrasts,xtxitx,svi,varargin)
% function [eff,eff_vector] = calcEfficiency(contrastweights,contrasts,xtxitx,svi,[D-opt flag])
%
% contrastweights   a row vector of length equal to number of rows of contrasts indicating relative weights
% contrasts         a matrix of contrasts to test, one row per contrast, one col per condition in X
% xtxitx            pinv(X), where X is the design matrix
% svi               S * Vi * S', where S is the filtering matrix and Vi is the intrinsic autocorrelation matrix
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
% by Tor Wager
% Documentation update: 9/18/09
% 
% Note: efficiency depends on scaling of contrastweights, contrasts, and
% hrf.  Comparing different models (e.g., with different basis sets) may be very tricky.
% in this implementation, contrastweights should sum to
% length(contrastweights) and be all positive for comparing designs.
%
% Power is proportional to sqrt(efficiency), holding constant the:
%  1) "optimal" effect size, magnitude over sigma
%  2) group sample size, and assuming between-subjects error = 0
% t = (M/sigma)/(1/sqrt(var)) =  (M/sigma) * sqrt(eff)
% assume M/sigma = 1, t = sqrt(eff)
% follows non-central t distribution with mean M/sigma
% p = 1 - tcdf(Z, df)
% p = 1 - tcdf(sqrt(eff), size(X, 1) - size(X, 2))
% -log(p) is -log of expected p-value given design and std. effect size

if isempty(svi), svi = eye(size(xtxitx,2)); end

Dflag = 0;
if length(varargin) > 0, Dflag = varargin{1}; end

if ~isempty(contrasts)											% compute the variance of the contrasts
    try
        eff_vector = diag(contrasts*xtxitx*svi*xtxitx'*contrasts');
     
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
            eff = 1 ./ inv(det(xtxitx)) .^ (1./size(xtxitx,2));   % HERE xtxitx SHOULD BE THE FISCHER INFO MTX X'X
        else
            eff = length(eff_vector)./(contrastweights * eff_vector);
        end
        eff_vector = 1 ./ eff_vector;
        % it's 1/ because we want to minimize the variance, but we maximize
        % fitnessMatrix value.

    catch
        disp('Error calculating efficiency: contrast weights and number of contrasts do not match!!!')
        disp(' I am trying to do this calculation:')
        disp('eff = length(eff_vector)./(contrastweights * eff_vector);')

        cwsize = size(contrastweights);
        effsize = size(eff_vector);

        disp(' ');
        disp('I am trying to multiply these:')
        fprintf('contrastweights (cwsize) * eff_vector (%01d x %01d) \n ', cwsize, effsize)
        disp('The most common causes are too many or too few entries in contrastweights.');
        error('Exiting.');
    end

else
    if isempty(contrastweights),
        contrastweights = ones(1,size(xtxitx,1)-1); % one for each beta, without intercept
    end

    try
        eff_vector = diag(xtxitx*svi*xtxitx');
        eff_vector = 1 ./ eff_vector;
        
        if Dflag
            eff = 1 ./ inv(det(xtxitx)) .^ (1./size(xtxitx,2));   % HERE xtxitx SHOULD BE THE FISCHER INFO MTX X'X
        else
            % variance of chosen regressors.
            eff = length(eff_vector)./([contrastweights 0] * eff_vector); 
            % eff_vector(eff_vector == 0) = .0001;
        end
        
    catch
        disp(['no contrasts version: wrong sizes!!!'])
        contrastweights = [contrastweights 0]
        whos xtxitx
        whos svi
        diag(xtxitx*svi*xtxitx')
        error('Exiting...')
    end

end

return