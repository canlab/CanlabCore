function [params_obj, hrf_obj, params_obj_dat, hrf_obj_dat, info] = hrf_fit(SPM,T,method,mode)
% HRF estimation on fmri_data class object
%
% HRF estimation function for a single voxel;
%
% Implemented methods include: IL-model (Deterministic/Stochastic), FIR
% (Regular/Smooth), and HRF (Canonical/+ temporal/+ temporal & dispersion).
% With SPM.mat, TR, and experimental design are imported.
%
% :Inputs:
%
%   **SPM**
%        SPM.mat file
%
%   **T**
%        length of estimated HRF ij seconds
%
%   **type**
%        Model type: 'FIR', 'IL', or 'CHRF'
%
%   **mode**
%        Mode
%
% :Model Types:
%
% A. **Fit HRF using IL-function**
%        Choose mode (deterministic/stochastic)
%           - 0 - deterministic aproach 
%           - 1 - simulated annealing approach
%
%        Please note that when using simulated annealing approach you
%        may need to perform some tuning before use.
%
% B. **Fit HRF using FIR-model**
%        Choose mode (FIR/sFIR)
%           - 0 - FIR
%           - 1 - smooth FIR
%
% C. **Fit Canonical HRF**
%        Choose mode (FIR/sFIR)
%           - 0 - FIR 
%           - 1 - smooth FIR        
%
%
% ..
%    Created by Michael Sun on 02/20/24
% ..

if isstring(SPM) || ischar(SPM)

    load(SPM);
    if ~exist('SPM', 'var')
        error('Passed in filepath is not an SPM.mat file.')
    end

end


if isstruct(SPM)
    fnames=unique({SPM.xY.VY.fname})';

    parfor i=1:numel(fnames)
        
        if contains(fnames{i}, '/') && ispc
            % PC Path conversion from Unix
            if strcmp(fnames{i}(1,1:2), '\\')
                d{i}=fmri_data(strrep(fnames{i}, '/', '\'));
            else
                d{i}=fmri_data(['\',strrep(fnames{i}, '/', '\')]);
            end

        else
            d{i}=fmri_data(fnames{i});
        end
    end

    % Transform timeseries data into the way SPM desires.
    gkwy_d=spmify(d,SPM);

    % Extract TR
    TR = SPM.xY.RT;

    if isempty(T)
        disp('T is empty, generating T as 2 times the maximum duration for each task regressor.');
        for i = 1:numel(SPM.Sess)
            T{i}=cellfun(@(cellArray) 2*ceil(max(cellArray)), {SPM.Sess(i).U.dur});
        end
    elseif numel(SPM.Sess)>1 && ~iscell(T)
        % error('SPM structures concatenating multiple runs must have time windows passed in with a matching number of cell arrays.');
        T=repmat({T}, 1, numel(SPM.Sess));
    elseif iscell(T) && numel(SPM.Sess)~=numel(T)
        error('Incompatible number of runs in SPM and in cell-array T.')
    end

    % %% ONE WAY
    % % run hrf_fit separately for every d
    for i=1:numel(d)
        % Reassign data.
        d{i}.dat=gkwy_d{i};

        % Extract TR
        TR = SPM.xY.RT;
        Runc=generateConditionTS(numel(SPM.Sess(i).row), [SPM.Sess(i).U.name], {SPM.Sess(i).U.ons}, {SPM.Sess(i).U.dur});
        
        % There's a need to pull apart the SPM design matrix for every run
        % Column of regressors + covariates for each session + intercept
        % X=SPM.xX.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i)]);

        % Trouble is, X is not filtered
        % X=spm_filter(SPM.xX.K(i), SPM.xX.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i)]));

        % Trouble here is, task regressors may have to be re-generated for
        % FIR and sFIR if the original design matrix is cHRF, so check the
        % basis function
        if strcmpi(method, 'FIR') && ~contains(SPM.xBF.name, 'FIR')
            disp('Passed in SPM structure is not FIR. Reconstructing task-regressors for FIR fit');
            
            len = numel(SPM.Sess(i).row);

            % Task regressors for each run can be found here:
            numstim=numel(SPM.Sess(i).U);

            % Make the design matrix:
            DX_all = cell(1, numstim); % Store DX matrices for each condition
            tlen_all = zeros(1, numstim); % Store tlen for each condition
            
            for s=1:numstim
                t = 1:TR:T{i}(s);
                tlen_all(s) = length(t);
                DX_all{s} = tor_make_deconv_mtx3(Runc(:,s), tlen_all(s), 1);
            end
            DX = horzcat(DX_all{:});
            
            % due to horzcat, we will have multiple intercepts in this design matrix
            % therefore, we'll find the intercepts, drop them, and add an intercept at the end of the matrix
            intercept_idx = find(sum(DX)==len);
            copyDX = DX;
            copyDX(:,intercept_idx) = [];
            DX = [copyDX ones(len,1)];
       
            % Covariate Design Matrix for each session (without intercept):
            NX=[SPM.Sess(i).C.C];

            % Concatenate the Task regressors with Covariates
            X=[DX, NX];
            % Filter
            X=spm_filter(SPM.xX.K(i), X);

        else
            % Otherwise set X to be the filtered design matrix of that
            % section.
            X=spm_filter(SPM.xX.K(i), SPM.xX.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i)]));
        end

        [params_obj{i}, hrf_obj{i}, params_obj_dat{i}, hrf_obj_dat{i}] = hrf_fit(d{i},TR,Runc,T{i},method,mode,X);

        info{i}.X=X;
        info{i}.names=[SPM.Sess(i).U.name]';

    end

end

end