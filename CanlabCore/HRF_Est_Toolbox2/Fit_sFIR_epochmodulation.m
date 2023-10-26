function [hrf, fit, e, param] = Fit_sFIR_epochmodulation(tc, TR, Run, T, mode)
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
% [10/26/2023 - HJ v1.2]: Updated the code to Allow for Different Epoch Lengths
% - handles varying epoch lengths with the addition of T vector, instead of a scalar T
% - made adjustments accordingly to handle T vector in variable `tlen_all`
% - based on varying epoch lengths, dynamically calculates matrix using `tor_make_deconv_mtx3.m`
% - (dynamic matrix is concatenated into one design matrix at the end)
% - made adjustments accordingly to the estimated hrf variable: instead of being a vector, hrf is now a matrix with varying lengths of hrf estimation, due to varying epoch length inputted via T.
% - lastly, varying lengths of hrf is concatenated into one matrix (double: longest length hrf estimation x number of conditions)

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

if mode == 1
   
    MRI = zeros(sum(tlen_all)+numstim); % adjust size based on varying tlen
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
    
    b = inv(DX'*DX+sig^2*MRI)*DX'*tc;
    fit = DX*b;
    e = tc - DX*b; 

elseif mode == 0

    b = pinv(DX)*tc;
    fit = DX*b;
    e = tc - DX*b;
    
end

numstim = length(tlen_all);
hrf = cell(1, numstim);

for i = 1:numstim
    hrf{i} = zeros(tlen_all(i), numstim);
end

%     hrf =zeros(tlen,numstim);
param = zeros(3,numstim);

for i=1:numstim
    hrf{:,i} = b(((i-1)*tlen_all(i)+1):(i*tlen_all(i)))';
    param(:,i) = get_parameters2(hrf{:,i},(1:tlen_all(i)));
end;

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
