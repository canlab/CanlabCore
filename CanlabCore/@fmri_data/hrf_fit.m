function [params_obj, hrf_obj, params_obj_dat, hrf_obj_dat] = hrf_fit(obj,TR,Runc,T,method,mode, varargin)
% HRF estimation on fmri_data class object
%
% HRF estimation function for a single voxel;
%
% Implemented methods include: IL-model (Deterministic/Stochastic), FIR
% (Regular/Smooth), and HRF (Canonical/+ temporal/+ temporal & dispersion)
%
% :Inputs:
%
%   **obj**
%        fMRI object or cell-array of fMRI-objects
%
%   **TR**
%        time resolution
%
%   **Runs**
%        expermental design
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
% :Examples:
% SIMULATE DATA AND RUN
% ::
%
%     %params for sim and fitting
%     TR = 2;   % repetition time (sec)
%     n = 200;  % time points measured (for simulation) must be multiple of 10
%     T = 30;   % duration of HRF to estimate (seconds)
%     nconds = 2; % num conditions
%     nevents = 8; % events per condition
%
%     % Create fake data
%     h = spm_hrf(TR);
%     y = zeros(n, 1);
%
%     % onsets - indicator
%     Condition = {};
%     for i = 1:nconds
%         Condition{i} = zeros(n,1);
%         wh = randperm(n);
%         Condition{i}(wh(1:nevents)) = 1;
%
%         ytmp{i} =  conv(Condition{i}, h);
%         ytmp{i} = ytmp{i}(1:n);
%     end
%
%     y = sum(cat(2, ytmp{:}), 2);
%
%     dat = fmri_data('VMPFC_mask_neurosynth.img');  % AVAILABLE ON WIKI IN MASK GALLERY
%     dat = threshold(dat, [5 Inf], 'raw-between');
%
%     v = size(dat.dat, 1); % voxels in mask
%     dat.dat = repmat(y',v, 1) + .1 * randn(v, n);
%
%     % Fit data - estimate HRFs across the brain mask
%     [params_obj hrf_obj] = hrf_fit(dat,TR, Condition, T,'FIR', 1);
%
%     hrf = fmri_data('HRF_timecourse_cond0001.img');
%     hrf = remove_empty(hrf);
%     create_figure('hrfs', 1, 2); 
%     plot(hrf.dat');
%     title('Condition 1')
%     hrf = fmri_data('HRF_timecourse_cond0002.img');
%     hrf = remove_empty(hrf);
%     subplot(1, 2, 2);
%     plot(hrf.dat');
%     title('Condition 2')
%
% ..
%    Created by Martin Lindquist on 04/11/14
% ..

% Added varargin to pass in a custom design matrix. - Michael Sun, Ph.D.
% 02/09/2024

fprintf('HRF estimation\n')

if iscell(obj)

    for i = numel(obj)
        
        
        [params_obj{i}, hrf_obj{i}, params_obj_dat{i}, hrf_obj_dat{i}] = hrf_fit(obj{i},TR,Runc,T,method,mode, varargin);


    end

end





% SPM=[];
if ~isempty(varargin)
    fprintf('Using passed in design matrix\n')
    % SPM = varargin{1};
    % 
    % if ischar(SPM) || isstr(SPM)
    %     load(SPM);
    % end
    % 
    % scan_nums=numel(SPM.Sess);

    % If you've got concatenated runs, you'll have to find a way to
    % identify which run your object is located in your SPM design matrix.
    % dat=spmify(obj, SPM, scan_num);
    % obj.dat=dat;
    DX=varargin{1};
    % task_regressors=varargin{2}; % Identify the columns of your task-regressors in DX.

end

%fhan = @(tc) hrf_fit_one_voxel(tc,TR,Runc,T,method,mode);

if ~isempty(varargin)
    fhan = @(tc) hrf_fit_one_voxel_meval(tc,TR,Runc,T,method,mode, varargin{:});
else
    [~,~,~,~,info]=hrf_fit_one_voxel(obj.dat(1,:)',TR,Runc,T,method,mode);
    fhan = @(tc) hrf_fit_one_voxel_meval(tc,TR,Runc,T,method,mode, info.DX, 'invertedDX', info.PX);
end

if isstruct(obj)
    SPM=obj;
    fnames=unique({SPM.xY.VY.fname})';

    for i=1:numel(fnames)
        
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


    hrf_fit(obj)


end




narg_to_request = length(Runc) * 2;
myoutargs = cell(1, narg_to_request);

[myoutargs{:}] = matrix_eval_function(obj.dat',fhan);

%[hrf_vals_brain param_vals_brain] = matrix_eval_function(obj.dat',fhan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numstim = length(Runc);  % number of conditions

shell_obj = mean(obj);  % shell object
    
for i=1:numstim

    % hrf_obj_i = shell_obj;
    % params_obj_i = shell_obj;

    hrf_obj_dat{i} = shell_obj;
    params_obj_dat{i} = shell_obj;

    % Generate fmri_dat objects
    params_obj{i}.dat = myoutargs{numstim + i};    
    hrf_obj_dat{i}.dat= myoutargs{i};

    %%%%%%%%%%%%%%%%%%%
    % Assign param data
    % params_obj_i.dat = myoutargs{numstim + i};
    
    %%%%%%%%%%%%%%%%%%%%
    % Assign HRF data
    % hrf_obj_i.dat = myoutargs{i};

    % params_obj_i.fullpath = sprintf('HRF_params_cond%04d.img', i); 
    % hrf_obj_i.fullpath = sprintf('HRF_timecourse_cond%04d.img', i);
    
    %%%%%%%%%%%%%%%%%%%%
    % Write data
    % write(params_obj_i);
    % write(hrf_obj_i);

end

hrf_obj = myoutargs(1:numstim);
params_obj = myoutargs(numstim+1 : end); % Estimated amplitude, height, and width of each condition


end


function varargout = hrf_fit_one_voxel_meval(tc,TR,Runc,T,method,mode, varargin)
% parse outputs into separate row vectors (for matrix_eval_function)
% outputs must be row vectors

if isempty(varargin)
    [~, ~, ~, ~, info] = hrf_fit_one_voxel(tc,TR,Runc,T,method,mode);
    [h, fit, e, param] = hrf_fit_one_voxel_meval(tc,TR,Runc,T,method,mode, info.DX, 'invertedDX', info.PX);
else
    % Use a passed-in custom design-matrix - Michael Sun, Ph.D. 02/13/2024
    [h, fit, e, param] = hrf_fit_one_voxel(tc,TR,Runc,T,method,mode, varargin{:});
end    

varargout = {};
for i = 1:size(h, 2)
    varargout{end+1} = h(:, i)';
end

for i = 1:size(param, 2) % 
    varargout{end+1} = param(:, i)';
end


end

