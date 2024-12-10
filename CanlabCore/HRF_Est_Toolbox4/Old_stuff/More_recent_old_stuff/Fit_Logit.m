function [hrf, fit, e, param] = Fit_Logit(tc,Run,t,mode)
% function [hrf, fit, e, param] = Fit_Logit(tc,Run,t,mode)
%
% Fits FIR and smooth FIR model  
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% mode  - deterministic or stochastic
%   options:
%       0 - deterministic aproach 
%       1 - simulated annealing approach
%       Please note that when using simulated annealing approach you
%       may need to perform some tuning before use.
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% Created by Martin Lindquist on 10/02/09
% Last edited: 05/26/10 (ML)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the Logit model

numstim = length(Run);
len = length(Run{1});

V0 = [ 1 6 1 0.5 10 1 15];                % initial values for logit fit
V0 = repmat(V0,1,numstim);

if (mode == 1 && numstim>1),
    disp('Multi-stimulus annealing function currently not implemented. Switching to "deterministic mode"')
    mode = 0;
end;

% Estimate theta (logit parameters)

if (mode == 1)
    disp('Stochastic Mode');
    
    Runs = Run{1};
    [theta,HH,C,P,hrf,fit,e,param] = Anneal_Logit(V0,t,tc,Runs);       
    
elseif (mode == 0)
    disp('Deterministic Mode');
    [theta, hrf, fit, e, param] = Det_Logit(V0,t,tc,Run);

end

end
