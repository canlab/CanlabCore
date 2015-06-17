function [eff,eff_vector] = calcEfficiency(contrastweights,contrasts,xtxitx,svi,varargin)
% function [eff,eff_vector] = calcEfficiency(contrastweights,contrasts,xtxitx,svi,[D-opt flag])
%
% contrastweights   a row vector of length equal to number of rows of contrasts indicating relative weights
% contrasts         a matrix of contrasts to test, one row per contrast, one col per condition in X
% xtxitx            pinv(X), where X is the design matrix
% svi               S * Vi, where S is the filtering matrix and Vi is the intrinsic autocorrelation matrix
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
%                   example: eff = calcEfficiency([],[],X'*X,[],1)  % d-optimality measure
%
% Another example: Try on autocorrelated and uncorrelated data
% X is a Fisher Info Matrix, XtX
% Vi = getV('make',[1 .5 .2 .1 0],10);
% X = eye(10)*2; eff = calcEfficiency([],[],X,Vi,1)
% X = eye(10)*2; eff = calcEfficiency([],[],X,[],1)
%
% by Tor Wager

if isempty(svi), svi = eye(size(xtxitx,2));,end

Dflag = 0;
if length(varargin) > 0, Dflag = varargin{1};, end

if ~isempty(contrasts)											% compute the variance of the contrasts
           try
              vcvm = contrasts*xtxitx*svi*xtxitx'*contrasts';   % var-cov matrix of contrasts
			  if Dflag
                  eff = det(vcvm) .^ (1./size(vcvm,2));
                  %same as 1 ./ inv(det(xtxitx)) .^ (1./size(xtxitx,2));   % HERE xtxitx SHOULD BE THE FISCHER INFO MTX X'X
              else
                  eff_vector = diag(vcvm);
                  eff = length(eff_vector)./(contrastweights * eff_vector);
                  eff_vector = 1 ./ eff_vector;
              end
          
           catch
             disp('using contrasts: matrices or contrasts are wrong sizes!!!')
             disp('Trying to do: eff_vector = diag(contrasts*xtxitx*svi*xtxitx''*contrasts'')')
	        disp('...and then eff = length(eff_vector)./(contrastweights * eff_vector);')

             contrasts = contrasts	
             whos xtxitx
             whos svi
	         contrastweights = contrastweights
             xtxitx = diag(xtxitx*svi*xtxitx')
             error('Exiting...')
           end
           
		   
 else	
           if isempty(contrastweights),
               contrastweights = ones(1,size(xtxitx,1)-1);
           end

           try
			   vcvm = xtxitx*svi*xtxitx';
               if Dflag
                  eff = det(vcvm) .^ (1./size(vcvm,2));   % HERE xtxitx SHOULD BE THE FISCHER INFO MTX X'X
                else
                    eff_vector = diag(vcvm);
                    eff = length(eff_vector)./([contrastweights 0] * eff_vector); % variance of chosen regressors.
                    eff_vector = 1 ./ eff_vector;
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