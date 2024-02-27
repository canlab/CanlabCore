function [hrf, fit, e, param, DXinfo] = Fit_sFIR_epochmodulation(tc, TR, Run, T, mode, varargin)
%
% Fits FIR and smooth FIR model  
%
% INPUTS:
% ----------
%   tc : (scan_length x 1) double vector 
%      time course of a single voxel.
%
%   TR : float
%       time resolution of one scan. 
%       e.g. use `0.5` for a scan collected every 0.5 seconds.
%
%   Runs : cell array (1 x num_condition)
%      experimental design matrix. each cell, the size of (1 x num_condition) contains 
%      a (scan_length x 1) double vector representing the design for that condition.
%
%   T : double array
%      desired length of the estimated HRF. indicates the number of beta 
%      coefficients to estimate, for each condition.
%      e.g. [30, 30, 30, 30, 30, 30, 20, 12] for a design with 8 conditions.
%
%   mode : int {0, 1}
%      options for FIR.
%      0 - Standard FIR.
%      1 - Smooth FIR.
%
% 
% OUTPUTS:
% -------
%   hrf : double array (time_point x num_condition)
%      estimated hemodynamic response function.
%      estimated HRF length differs depending on inputted T, across conditions.
%      to account for this difference, the hrf matrix is the size of the largest time 
%      points estimated, which is max(T/TR). the shorter events are filled with NaNs
%
%   fit : double array (scan_length x 1) 
%      estimated time course.
%
%   e : double array (scan_length x 1) 
%      residual time course.
%
%   param : double array (3 x num_condition)
%      estimated amplitude, height and width
%
% Created by Martin Lindquist on 10/02/09
% Last edited: 05/17/13 (ML)
% 
% [10/27/2023 - HJ v1.2]: BUGFIX: remove redundant intercept due to multiple design matrix estimation
% - due to horzcat, we will have multiple intercepts in this design matrix
% - therefore, we find the intercepts, drop them, and add an intercept at the end of the matrix

% [10/26/2023 - HJ v1.2]: Updated the code to Allow for Different Epoch Lengths
% - handles varying epoch lengths with the addition of T vector, instead of a scalar T
% - made adjustments accordingly to handle T vector in variable `tlen_all`
% - based on varying epoch lengths, dynamically calculates matrix using `tor_make_deconv_mtx3.m`
% - (dynamic matrix is concatenated into one design matrix at the end)
% - made adjustments accordingly to the estimated hrf variable: instead of being a vector, hrf is now a matrix with varying lengths of hrf estimation, due to varying epoch length inputted via T.
% - lastly, varying lengths of hrf is concatenated into one matrix (double: longest length hrf estimation x number of conditions)

% [2/19/2024 - MS v2.0]: Updated the code to output design-matrix information, and allow input of multiple tc,
% pass in a custom design matrix, and to pass in an SPM structure

% Debugging Hack, import your own Design Matrix from SPM - MS
if ~isempty(varargin)
    if isstruct(varargin{1})
        % 1. Decide what are regressors and what are covariates here
        % This is very difficult since I will have to pass in g, K, and W
        % all from SPM.mat
        
        % So ultimately we will want to concatenate the FIR task regressors
        % and the unfiltered covariates, and then filter the whole thing.
        
        SPM=varargin{1};

        if numel(SPM.Sess)>1 && ~iscell(tc) % If SPM has concatenated runs but there's only one time series
            error('SPM structures concatenating multiple runs must have timecourses passed in with a matching number of cell arrays');
        elseif iscell(tc) && numel(SPM.Sess)~=numel(tc)
            error('Incompatible number of runs in SPM and in cell-array tc')
        end

        if isempty(T)
            for i = 1:numel(SPM.Sess)
                T{i}=cellfun(@(cellArray) 2*ceil(max(cellArray)), {SPM.Sess(i).U.dur});
            end
        elseif numel(SPM.Sess)>1 && ~iscell(T)
            error('SPM structures concatenating multiple runs must have time windows passed in with a matching number of cell arrays.');
        elseif iscell(T) && numel(SPM.Sess)~=numel(T)
            error('Incompatible number of runs in SPM and in cell-array T.')
        end

        for s = 1:numel(SPM.Sess)
            len = numel(SPM.Sess(s).row);

            % Task regressors for each run can be found here:
            numstim=numel(SPM.Sess(s).U);
            TR = SPM.xY.RT;

            % Make the design matrix:
            DX_all = cell(1, numstim); % Store DX matrices for each condition
            tlen_all = zeros(1, numstim); % Store tlen for each condition
            Runs{s}=generateConditionTS(numel(SPM.Sess(s).row), [SPM.Sess(s).U.name], {SPM.Sess(s).U.ons}, {SPM.Sess(s).U.dur});
            
            for i=1:numstim
                t = 1:TR:T{s}(i);
                tlen_all(i) = length(t);
                DX_all{i} = tor_make_deconv_mtx3(Runs{s}(:,i), tlen_all(i), 1);
            end
            DX = horzcat(DX_all{:});
            
            % due to horzcat, we will have multiple intercepts in this design matrix
            % therefore, we'll find the intercepts, drop them, and add an intercept at the end of the matrix
            intercept_idx = find(sum(DX)==len);
            copyDX = DX;
            copyDX(:,intercept_idx) = [];
            DX = [copyDX ones(len,1)];
       
            % Covariate Design Matrix for each session (without intercept):
            NX=[SPM.Sess(s).C.C];

            % Concatenate the Task regressors with Covariates
            X=[DX, NX];
            % Filter
            DX_cov{s}=spm_filter(SPM.xX.K(s), X);
        end

        if numel(DX_cov) > 1
            for i=1:numel(DX_cov)
                [hrf{i}, fit{i}, e{i}, param{i}, DXinfo{i}]=Fit_sFIR_epochmodulation(tc{i}, TR, Runs{i}, T{i}, mode, DX_cov{i});
            end
            return;
        end

    else
        % If not an SPM struct, allow a design matrix to be passed in.
        DX_cov=varargin{1};

        numstim = length(Run);
        tlen_all = zeros(1, numstim); % Store tlen for each condition
        for i=1:numstim
            t = 1:TR:T(i);
            tlen_all(i) = length(t);
        end
    end
    

    if mode == 1
        % pen = {toeplitz([1 .3 .1])}
        % n_nuis = 20;  % num of nuisance
        % pen{2} = zeros(20); % for nuisance, add no penalty
        % blkdiag(pen{:})

        MRI = zeros(sum(tlen_all)+1); % adjust size based on varying tlen and intercept at end
        start_idx = 1;
    
        for i=1:numstim
            tlen = tlen_all(i); % get the tlen for this stimulus
    
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
        cov_num = size(DX_cov,2)-size(MRI,2);
        % pen = sig^2*MRI; % Regularization Penalty Matrix
        % pen=pad(pad(pen, cov_num)',cov_num)';
        % pen{2} = zeros(20); % for nuisance, add no penalty
        pen{1} = sig^2*MRI; % Regularization Penalty Matrix
        pen{2} = zeros(cov_num); % for nuisance, add no penalty
        pen=blkdiag(pen{:});

        % disp(pen) % Check what this looks like.
        % X = [DX, NX]; % Full Design Matrix

        b = inv(DX_cov'*DX_cov+pen)*DX_cov'*tc;
        % b = DX'*tc;
        fit = DX_cov*b;
        e = tc - DX_cov*b; 
        DXinfo.DX=DX_cov';

    elseif mode == 0
        b = pinv(DX_cov)*tc;
        fit = DX_cov*b;
        e = tc - DX_cov*b;
        DXinfo.DX=DX_cov;

    end

else

    % If nothing is passed in varargin:

    
    numstim = length(Run);
    len = length(Run{1});
    
    Runs = zeros(len,numstim);
    for i=1:numstim
        Runs(:,i) = Run{i};
    end
    
    DX_all = cell(1, numstim); % Store DX matrices for each condition
    tlen_all = zeros(1, numstim); % Store tlen for each condition
    
    
    for i=1:numstim
        t = 1:TR:T(i);
        tlen_all(i) = length(t);
        DX_all{i} = tor_make_deconv_mtx3(Runs(:,i), tlen_all(i), 1);
    end
    DX = horzcat(DX_all{:});
    
    % due to horzcat, we will have multiple intercepts in this design matrix
    % therefore, we'll find the intercepts, drop them, and add an intercept at the end of the matrix
    intercept_idx = find(sum(DX)==len);
    copyDX = DX;
    copyDX(:,intercept_idx) = [];
    DX = [copyDX ones(len,1)];
    
    if mode == 1
    
        MRI = zeros(sum(tlen_all)+1); % adjust size based on varying tlen and intercept at end
        start_idx = 1;
    
        for i=1:numstim
            tlen = tlen_all(i); % get the tlen for this stimulus
    
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
    
        % DX
        % 
        % NX = ;% Nuisance covariates and intercepts
    
        % X = [DX, NX]; % Full Design Matrix
        pen = sig^2*MRI; % Regularization Penalty Matrix
        % disp(pen) % Check what this looks like.
        cov_num = size(DX,2)-size(pen, 2);
    
        pen=pad(pad(pen, cov_num)',cov_num)';
    
        % multiply a 0 penalty with NX.
    
        % gKW the DX Design Matrix
    
    
        % b = inv(DX'*DX+sig^2*MRI)*DX'*tc;
        b = inv(DX'*DX+pen)*DX'*tc;
        fit = DX*b;
        e = tc - DX*b; 
    
        % Code added by Michael Sun, PhD 02/05/2024
        DXinfo.DX=inv(DX'*DX+sig^2*MRI)*DX'; % This DX is essentially regression coefficients
        DXinfo.DX=DXinfo.DX';
    
        % DXinfo.DX=DX;
        DXinfo.sig=sig;
        DXinfo.MRI=MRI;
    
    elseif mode == 0
    
        b = pinv(DX)*tc;
        fit = DX*b;
        e = tc - DX*b;
    
        % Coded added by Michael Sun, Ph.D.
        DXinfo.DX=DX;
    
    end



end

numstim = length(tlen_all);
hrf = cell(1, numstim);

for i = 1:numstim
    hrf{i} = zeros(tlen_all(i), numstim);
end

%     hrf =zeros(tlen,numstim);
param = zeros(3,numstim);

for i=1:numstim
    % hrf{:,i} = b(((i-1)*tlen_all(i)+1):(i*tlen_all(i)))'; % Buggy line Edit is below: Michael Sun, 10/27/2023   
    start_idx = sum(tlen_all(1:i-1)) + 1;
    end_idx = start_idx + tlen_all(i) - 1;
    hrf{:,i} = b(start_idx:end_idx)';
    param(:,i) = get_parameters2(hrf{:,i}, (1:tlen_all(i)));

    DXinfo.tlen{i}=[start_idx:end_idx];
end

% HJ: concatenate estimated hrf
% since different event lengths have different lengths of hrf estimation, 
% we'll identify the maximum length of hrf estimation timepoints; 
% for the shorter events, we'll pad them with zeros.
maxLen = max(cellfun(@length, hrf));
resultMatrix = NaN(maxLen, numstim); % Pre-allocate a matrix of NaN values
for i = 1:numel(hrf)
    resultMatrix(1:length(hrf{i}), i) = hrf{i}';
end
hrf = resultMatrix;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  validation plots: design matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% imagesc(DX); % Display matrix as an image
% colorbar; % Add colorbar to the figure
% title('Heatmap of DX');
% xlabel('Columns of DX');
% ylabel('Rows of DX');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % validation estimated curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numCells = length(hrf)';
% % Initialize a figure
% figure;
% for i = 1:numCells
%     % Create a subplot for each cell
%     subplot(2, 4, i);  % Assuming a 2x4 grid for 8 panels
    
%     % Plot the values of the double array
%     plot(hrf{i}, 'o-');
    
%     % Label axes and provide title for each subplot
%     xlabel('# of time points');
%     ylabel('Voxel-wise Amplitude');
%     title([num2str(length(hrf{i})) ' time points estimated']);  % Updated title
    
%     % Adjust y-axis for better visualization
%     ylim([min(hrf{i})-0.1, max(hrf{i})+0.1]);
% end
