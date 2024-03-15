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
        
        if ~contains(SPM.xBF.name, 'Finite Impulse Response')
            % For now, assume 'Canonical HRF'
            disp('T is empty, generating T as 2 times the maximum duration for each task regressor.');
            for i = 1:numel(SPM.Sess)
                T{i}=cellfun(@(cellArray) 2*ceil(max(cellArray)), {SPM.Sess(i).U.dur});
            end
        else
            
            for i = 1:numel(SPM.Sess)
                % T{i}=[SPM.xBF.order]; 
                % In normal SPM, the time window (length) and resolution (order) is set the same for all regressors.
                % T should be in seconds, it gets converted to TRs later.
                T{i}=repmat(SPM.xBF.length, size({SPM.Sess(i).U.dur}));

                % T{i}=repmat(SPM.xBF.order*TR, size({SPM.Sess(i).U.dur}));

            end
        end
    elseif numel(SPM.Sess)>1 && ~iscell(T)
        % error('SPM structures concatenating multiple runs must have time windows passed in with a matching number of cell arrays.');
        T=repmat({T}, 1, numel(SPM.Sess));
    elseif iscell(T) && numel(SPM.Sess)~=numel(T)
        error('Incompatible number of runs in SPM and in cell-array T.')
    end

    % %% ONE WAY
    % % run hrf_fit separately for every d

    % Tor: Better to pass in pseudoinverse matrix to speed up computation:
    parfor i=1:numel(d)
        % Reassign data.
        d{i}.dat=gkwy_d{i};

        % Extract TR
        TR = SPM.xY.RT;
        % 
        Runc=generateConditionTS(numel(SPM.Sess(i).row), [SPM.Sess(i).U.name], {SPM.Sess(i).U.ons}, {SPM.Sess(i).U.dur}, TR);
        
        % There's a need to pull apart the SPM design matrix for every run
        % Column of regressors + covariates for each session + intercept
        % X=SPM.xX.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i)]);

        % Trouble is, X is not filtered
        % X=spm_filter(SPM.xX.K(i), SPM.xX.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i)]));

        % Trouble here is, task regressors may have to be re-generated for
        % FIR and sFIR if the original design matrix is cHRF, so check the
        % basis function
        if strcmpi(method, 'FIR') && ~contains(SPM.xBF.name, 'Finite Impulse Response')
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

            % Design Matrix Regressors + intercept:
            DX = [copyDX ones(len,1)];
       
            % Covariate Design Matrix for each session (without intercept):
            NX=[SPM.Sess(i).C.C];

            % Concatenate the Task regressors with Covariates
            X=[DX, NX];
            % Filter
            % X=spm_filter(SPM.xX.K(i), X);

            % filter is probably performed with spm_sp('Set', xX.K*xX.W*xX.X)
            % X=spm_sp('Set', SPM.xX.K(i)*SPM.xX.W*SPM.xX.X)

            % xX.KxXs        = spm_sp('Set',spm_filter(xX.K,W*xX.X));    % KWX

            % Ke says preprocess the nuisance regressors first before
            % fitting the FIR, otherwise the design matrix has problematic
            % VIFs between the FIR indicators and the spike regressors.
            % xKNXs = spm_sp('Set',spm_filter(SPM.xX.K(i), SPM.xX.W(SPM.Sess(i).row, SPM.Sess(i).row)*NX));    % KW*Nuisance
            % NX    = full(xKNXs.X);
            % d{i}=canlab_connectivity_preproc(d{i}, 'additional_nuisance', NX, TR,'no_plots')
            % preproc_d=canlab_connectivity_preproc(d{i}, 'additional_nuisance', NX)

            % preproc_d=d{i};
            % d{i}.covariates=NX;
            % d{i}=preprocess(d{i}, 'resid',1); % Residualize out noise regressors

            % sFIR
            if mode == 1
                for s=1:numstim
                    t = 1:TR:T{i}(s);
                    tlen_all(s) = length(t);
                end

                MRI = zeros(sum(tlen_all)+1); % adjust size based on varying tlen and intercept at end
                start_idx = 1;
            
                for s=1:numstim
                    tlen = tlen_all(s); % get the tlen for this stimulus
            
                    C = (1:tlen)'*(ones(1,tlen));
                    h = sqrt(1/(7/TR));
            
                    v = 0.1;
                    sig = 1;
            
                    R = v*exp(-h/2*(C-C').^2);
                    RI = inv(R);
            
                    % Adjust the indices to account for varying tlen
                    end_idx = start_idx + tlen - 1;
                    MRI(start_idx:end_idx, start_idx:end_idx) = RI;
            
                    start_idx = end_idx + 1; % update the starting index for next iteration
                end
        
                % multiply a 0 penalty with NX.
                cov_num = size(X,2)-size(MRI,2); % Number of nuisance covariates
                pen{1} = sig^2*MRI; % Regularization Penalty Matrix
                pen{2} = zeros(cov_num); % for nuisance, add no penalty
                pen=blkdiag(pen{:});
        
                % disp(pen) % Check what this looks like.

                % Filter everything, and then perform the pseudoinverse
                xKXs = spm_sp('Set',spm_filter(SPM.xX.K(i), SPM.xX.W(SPM.Sess(i).row, SPM.Sess(i).row)*X));    % KW*Design
                X    = full(xKXs.X);

                PX = inv(X'*X+pen)*X';
                % PX = spm_sp('x-', X'*X+pen); % SPM's way of performing
                % the pseudoinversion. Need a way to put in in SPM's xKXs
                % space structure.
                % PX = PX*X'
            
            
            else
                % Filter and invert
                xKXs = spm_sp('Set',spm_filter(SPM.xX.K(i), SPM.xX.W(SPM.Sess(i).row, SPM.Sess(i).row)*X));    % KW*Design
                X    = full(xKXs.X);
                % PX   = pinv(X);
                PX = spm_sp('x-', xKXs); % SPM's way of performing the pseudoinversion.
                % figure, imagesc(X), clim([-.01 .01])
            end

        elseif contains(method, 'CHRF') && ~contains(SPM.xBF.name, 'hrf')
                % Problem here is we have onsets but no durations from a
                % FIR matrix. If we use FIR

                % Generate a design matrix
                Run=Runc;
                d = length(Run);
                len = length(Run{1});
                t=1:TR:T;
                
                % Constructing the Design Matrix X:
                X = zeros(len,p*d);
                param = zeros(3,d);
                
                [h, dh, dh2] = CanonicalBasisSet(TR);
                
                for j=1:d,
                    v = conv(Run{j},h);
                    X(:,(j-1)*p+1) = v(1:len);
                
                    % Computing the first derivative
                    if strcmpi(method, 'CHRF1')
                        v = conv(Run{j},dh);
                        X(:,(j-1)*p+2) = v(1:len);
                    end
                
                    % Computing the second derivative
                    if strcmpi(method, 'CHRF2')
                        v = conv(Run{j},dh2);
                        X(:,(j-1)*p+3) = v(1:len);
                    end
                end
        else
            % Otherwise set X to be the filtered design matrix of that
            % section.

            % If using sFIR, pre-penalize the pseudoinverted design matrix
            if strcmpi(method, 'FIR') && mode==1

                % Initialize an array to hold the length of regressors for each task
                tlen_all = [];
                
                % Extract condition names for the session
                condition_names = {SPM.Sess(i).U.name};
                
                % Loop through condition names to find lengths for each task
                for c = 1:length(condition_names)
                    condition = condition_names{c}{1};  % Adjusted for direct use without {1}
                    % Find indices of regressors matching this condition
                    condition_idx = find(contains(SPM.xX.name, ['Sn(', num2str(i), ') ', condition]));
                    
                    if ~isempty(condition_idx)
                        % Calculate the length of this task's block of regressors
                        task_length = max(condition_idx) - min(condition_idx) + 1;
                        % Add this length to the tlen_all array
                        tlen_all(end+1) = task_length;
                    end
                end

                MRI = zeros(sum(tlen_all)+1); % adjust size based on varying tlen and intercept at end
                start_idx = 1;
            
                for s=1:numstim
                    tlen = tlen_all(s); % get the tlen for this stimulus
            
                    C = (1:tlen)'*(ones(1,tlen));
                    h = sqrt(1/(7/TR));
            
                    v = 0.1;
                    sig = 1;
            
                    R = v*exp(-h/2*(C-C').^2);
                    RI = inv(R);
            
                    % Adjust the indices to account for varying tlen
                    end_idx = start_idx + tlen - 1;
                    MRI(start_idx:end_idx, start_idx:end_idx) = RI;
            
                    start_idx = end_idx + 1; % update the starting index for next iteration
                end

                % Create Penalty Matrix
                cov_num = size(SPM.Sess(i).C.C, 2);
                pen{1} = sig^2*MRI; % Regularization Penalty Matrix
                pen{2} = zeros(cov_num); % for nuisance, add no penalty
                pen=blkdiag(pen{:});
                % disp(pen) % Check what this looks like.
        
                % Invert
                X=SPM.xX.xKXs.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i), SPM.xX.iG]);

                % Now penalize SPM's task regressors
                PX = inv(X'*X+pen)*X';

            else
                X=SPM.xX.xKXs.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i), SPM.xX.iG]);
                % Troubleshooting:
                % For some reason the above differs from this:
                % X=spm_filter(SPM.xX.K(i), SPM.xX.X(SPM.Sess(i).row, [SPM.Sess(i).col, SPM.xX.iB(i), SPM.xX.iG]));
                % X_SPM=SPM.xX.xKXs.X(SPM.Sess(i).row,  [SPM.Sess(i).col, SPM.xX.iB(i), SPM.xX.iG]);
                % figure, imagesc(X_SPM), clim([-.01 .01])
                % Grab SPM's prefiltered pseudoinverse matrix
                PX=SPM.xX.pKX([SPM.Sess(i).col, SPM.xX.iB(i), SPM.xX.iG], SPM.Sess(i).row);
            
            end                        


        end

        % ~ 40 min
        % [params_obj{i}, hrf_obj{i}, params_obj_dat{i}, hrf_obj_dat{i}] = hrf_fit(d{i},TR,Runc,T{i},method,mode,X);
        % The reason this takes so long is because of the pseudoinversion
        % of the design matrix. If we have access to the SPM design matrix,
        % we should pass in both DX and its inversion.

        [params_obj{i}, hrf_obj{i}, params_obj_dat{i}, hrf_obj_dat{i}] = hrf_fit(d{i},TR,Runc,T{i},method,mode,X, 'invertedDX', PX);


        % Test

        % [params_obj1{i}, hrf_obj1{i}, params_obj_dat1{i}, hrf_obj_dat1{i}] = hrf_fit(d{i},TR,Runc,T{i},method,mode,X);
        % 
        % [params_obj2{i}, hrf_obj2{i}, params_obj_dat2{i}, hrf_obj_dat2{i}] = hrf_fit(d{i},TR,Runc,T{i},method,mode,X_SPM);
        % 
        % [params_obj3{i}, hrf_obj3{i}, params_obj_dat3{i}, hrf_obj_dat3{i}] = hrf_fit(d{i},TR,Runc,T{i},method,mode,X_SPM);
        % 
        % [params_obj4{i}, hrf_obj4{i}, params_obj_dat4{i}, hrf_obj_dat4{i}] = hrf_fit(d{i},TR,Runc,T{i},method,1,X_SPM);

        info{i}.X=X;
        info{i}.names=[SPM.Sess(i).U.name]';

        % Testing
        % Assuming i is defined, and your variables d, TR, Runc, T, method, mode, X, and X_SPM are initialized

        % Start or get the current parallel pool
        % pool = gcp();
        % 
        % % Asynchronously call hrf_fit with the first set of arguments
        % future1 = parfeval(pool, @hrf_fit, 4, d{i}, TR, Runc, T{i}, method, mode, X);
        % 
        % % Asynchronously call hrf_fit with the second set of arguments
        % future2 = parfeval(pool, @hrf_fit, 4, d{i}, TR, Runc, T{i}, method, mode, X_SPM);
        % 
        % % Wait for the first job to finish and collect its results
        % [params_obj1{i}, hrf_obj1{i}, params_obj_dat1{i}, hrf_obj_dat1{i}] = fetchOutputs(future1);
        % 
        % % Wait for the second job to finish and collect its results
        % [params_obj2{i}, hrf_obj2{i}, params_obj_dat2{i}, hrf_obj_dat2{i}] = fetchOutputs(future2);

    end
    % Possibly need to drop null regressors? No


    % betas = filenames(fullfile(pwd, 'data', 'sub-SID000743', '*heat_start*', '*ses-12*bodymap_run*.nii'));
    % apply_nps([betas(3), betas(4), betas(2), betas(1)])
    % 
    % % Bug testing
    % nps_test_cHRF={apply_nps(hrf_obj_dat3{1}), apply_nps(hrf_obj_dat3{2}), apply_nps(hrf_obj_dat3{3}), apply_nps(hrf_obj_dat3{4})}
    % nps_test_FIR={apply_nps(hrf_obj_dat1{1}), apply_nps(hrf_obj_dat1{2}), apply_nps(hrf_obj_dat1{3}), apply_nps(hrf_obj_dat1{4})}
    % nps_test_sFIR={apply_nps(hrf_obj_dat2{1}), apply_nps(hrf_obj_dat2{2}), apply_nps(hrf_obj_dat2{3}), apply_nps(hrf_obj_dat2{4})}
    % 
    % figure
    % plot([nps_test_FIR{1}{4}], 'r--'), hold on
    % plot([nps_test_FIR{2}{4}], 'g--'), hold on
    % plot([nps_test_FIR{3}{4}], 'b--'), hold on
    % plot([nps_test_FIR{4}{1}], 'm--'), hold on
    % 
    % plot([nps_test_FIR{1}{4}, nps_test_FIR{2}{4},nps_test_FIR{3}{4},nps_test_FIR{4}{1}]), hold on
    % 
    % plot([nps_test_sFIR{1}{4}], 'r-'), hold on
    % plot([nps_test_sFIR{2}{4}], 'g-'), hold on
    % plot([nps_test_sFIR{3}{4}], 'b-'), hold on
    % plot([nps_test_sFIR{4}{1}], 'm-'), hold on
    % 
    % plot([nps_test_sFIR{1}{4}, nps_test_sFIR{2}{4},nps_test_sFIR{3}{4},nps_test_sFIR{4}{1}]), hold on
    % 
    % hline(0,'k-')
    % legend({'run1', 'run2', 'run3'})
    % legend({'run1', 'run2', 'run3', 'run4'})
    % 
    % % compare with betas
    % 
    % % Check out what d looks like
    % nps_d=apply_nps(d)
    % plot([nps_d{:}])
    % 
    % % Checkout what spmify did
    % figure
    % nps_gkwyd=apply_nps(gkwy_data)
    % plot([nps_gkwyd{:}])
    % 
    % % Plot to compare
    % plot([nps_d{:}]), hold on
    % plot([nps_gkwyd{:}]), hold on
    % 
    % plot([nps_d{1}]), hold on
    % plot([nps_gkwyd{1}]), hold on
    % legend({'pre', 'post'})
    % 
    % plot([nps_d{2}]), hold on
    % plot([nps_gkwyd{2}]), hold on
    % legend({'pre', 'post'})
    % 
    % plot([nps_d{3}]), hold on
    % plot([nps_gkwyd{3}]), hold on
    % legend({'pre', 'post'})
    % 
    % plot([nps_d{4}]), hold on
    % plot([nps_gkwyd{4}]), hold on
    % legend({'pre', 'post'})
    % 







end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, dh, dh2] = CanonicalBasisSet(TR)

len = round(30/TR);
xBF.dt = TR;
xBF.length= len;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);

v1 = xBF.bf(1:len,1);
v2 = xBF.bf(1:len,2);
v3 = xBF.bf(1:len,3);

h = v1;
dh =  v2 - (v2'*v1/norm(v1)^2).*v1;
dh2 =  v3 - (v3'*v1/norm(v1)^2).*v1 - (v3'*dh/norm(dh)^2).*dh;

h = h./max(h);
dh = dh./max(dh);
dh2 = dh2./max(dh2);

end