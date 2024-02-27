function [h, fit, e, param] =  hrf_fit_one_voxel(tc,TR,Runc,T,method,mode, varargin)
% function [h, fit, e, param] =  hrf_fit_one_voxel(tc,TR,Runc,T,type,mode)
%
% HRF estimation function for a single voxel;
%
% Implemented methods include: IL-model (Deterministic/Stochastic), FIR
% (Regular/Smooth), and HRF (Canonical/+ temporal/+ temporal & dispersion)
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% type  - Model type
% mode  - Mode
%
% Options: p=1 - only canonical HRF
%          p=2 - canonical + temporal derivative
%          p=3 - canonical + time and dispersion derivative
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
%
% Created by Martin Lindquist on 04/11/14
%
%   Added a varargin to accept a design matrix to modulate the model fit.
%   - Michael Sun, Ph.D. 02/09/2024


if ~isempty(varargin)
    % Load in SPM structure if passed in
    % timecourse tc must now be gKWY to match SPM output.
    % Filtered Design Matrix SPM.xX.xKXs.X must be passed in instead of
    % allowing Fit_* to generate its own design matrix.

    % if ischar(varargin{1}) || isstring(varargin{1})
    %     load(varargin{1})
    % elseif isstruct(varargin{1})
    % 
    %     SPM=varargin{1};
    % end
    % 
    % SPM.xX.xKXs.X

    SPM_DX=varargin{1};

end


if (strcmp(method,'IL')),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using IL-function

% Choose mode (deterministic/stochastic)

% mode = 0;     % 0 - deterministic aproach 
                % 1 - simulated annealing approach
                % Please note that when using simulated annealing approach you
                % may need to perform some tuning before use.

            if exists('SPM')
                disp('IL-function with custom Design Matrix not yet implemented')
                [h, fit, e, param] = Fit_Logit2(tc,TR,Runc,T,mode);
            else
                [h, fit, e, param] = Fit_Logit2(tc,TR,Runc,T,mode);
            end

            param(2:3,:) = param(2:3,:).*TR;


elseif (strcmp(method,'FIR')),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using FIR-model

% Choose mode (FIR/sFIR)

% mode = 1;     % 0 - FIR 
                % 1 - smooth FIR
            
            % [h, fit, e, param] = Fit_sFIR(tc,TR,Runc,T,mode);
            

            if exist('SPM_DX', 'var')
                [h, fit, e, param] = Fit_sFIR_epochmodulation(tc,TR,Runc,T,mode, SPM_DX);

            else
                [h, fit, e, param] = Fit_sFIR_epochmodulation(tc,TR,Runc,T,mode);
            end



            param(2:3,:) = param(2:3,:).*TR;


elseif (strcmp(method,'CHRF')),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

% mode = 0;     % 0 - HRF 
                % 1 - HRF + temporal derivative
                % 2 - HRF + temporal and dispersion derivative
            
            p = mode + 1;


            if exist('SPM_DX', 'var')
                [h, fit, e, param] = Fit_Canonical_HRF(tc,TR,Runc,T,p,SPM_DX);
            else
                [h, fit, e, param] = Fit_Canonical_HRF(tc,TR,Runc,T,p);
            end


            param(2:3,:) = param(2:3,:).*TR;

else
                warning('Incorrect Model Type. Use IL, FIR or CHRF');
end



end % function

