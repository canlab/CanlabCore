function [h, fit, e, param] =  hrf_fit_one_voxel(tc,TR,Runc,T,method,mode)
% HRF estimation function for a single voxel;
%
% :Usage:
% ::
%
%     function [h, fit, e, param] =  hrf_fit_one_voxel(tc,TR,Runc,T,type,mode)
%
% Implemented methods include: IL-model (Deterministic/Stochastic), FIR
% (Regular/Smooth), and HRF (Canonical/+ temporal/+ temporal & dispersion)
%
% :Inputs:
%
%   **tc:**
%        time course
%
%   **TR:**
%        time resolution
%
%   **Runs:**
%        expermental design
%
%   **T:**
%        length of estimated HRF
%
%   **type:**
%        Model type
%
%   **mode:**
%        Mode
%
% :Options:
%    - p=1 - only canonical HRF
%    - p=2 - canonical + temporal derivative
%    - p=3 - canonical + time and dispersion derivative
%
% :Outputs:
%
%   **hrf:**
%        estimated hemodynamic response function
%
%   **fit:**
%        estimated time course
%
%   **e:**
%        residual time course
%
%   **param:**
%        estimated amplitude, height and width
%
% ..
%    Created by Martin Lindquist on 04/11/14
% ..


if (strcmp(method,'IL')),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using IL-function

% Choose mode (deterministic/stochastic)

% mode = 0;     % 0 - deterministic aproach 
                % 1 - simulated annealing approach
                % Please note that when using simulated annealing approach you
                % may need to perform some tuning before use.

            [h, fit, e, param] = Fit_Logit2(tc,TR,Runc,T,mode);
            param(2:3,:) = param(2:3,:).*TR;


elseif (strcmp(method,'FIR')),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using FIR-model

% Choose mode (FIR/sFIR)

% mode = 1;     % 0 - FIR 
                % 1 - smooth FIR
            
            [h, fit, e, param] = Fit_sFIR(tc,TR,Runc,T,mode);
            param(2:3,:) = param(2:3,:).*TR;


elseif (strcmp(method,'CHRF')),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

% mode = 0;     % 0 - HRF 
                % 1 - HRF + temporal derivative
                % 2 - HRF + temporal and dispersion derivative
            
            p = mode + 1;            
            [h, fit, e, param] = Fit_Canonical_HRF(tc,TR,Runc,T,p);
            param(2:3,:) = param(2:3,:).*TR;

else
                warning('Incorrect Model Type. Use IL, FIR or CHRF');
end



end % function

