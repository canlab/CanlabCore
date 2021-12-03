function single_trial_estimates(obj, fmri_data_obj)
% Write single trial estimates associated with an estimated fmri_model object.
% must have estimated the model (robustfit(obj); see fmri_model.robustfit)
% and saved hrf*.img images for each condition.
%
% Also input an fmri_data object with time series data.
%
% This function writes images, one 4-D image for each condition, with the
% number of frames equalling the number of trials (onsets) for that
% condition.
%
% It does this by constructing a separate design matrix for each voxel,
% which is based on the HRF estimates for that voxel for each condition.
% Fits for all conditions are added to the same model, so that their
% colinearity influences the single-trial parameter estimates.

c = obj.xX.cond_assignments;
nconds = size(c, 2);


% check fmri_data_obj data
% ---------------------------------------------------------------------
Y = fmri_data_obj.X;

fmri_data_obj.X = []; % save space

fprintf('Checking data');
badvox = any(isnan(Y)) | all(Y == 0);
anybadvox = any(badvox); % use later as well
if anybadvox
    fprintf('Voxels with bad values: %3.0f. (all zeros or some NaNs)\n', sum(badvox));
    Y(:, badvox) = [];
    % we will have to insert these later.
end


% load hrf images
% ---------------------------------------------------------------------
disp('Loading mask.img from current (analysis) directory:');
mask = fmri_mask_image('mask.img');

disp('Single trial estimates for fmri_model object');
disp('Loading hrf estimates from images...');

for i = 1:nconds
    
    fprintf('Condition %3.0f, %s', i, obj.xX.name{i});
    
    hrfname = fullfile(pwd, sprintf('hrf_%s_cond%03d.img', obj.xX.name{i}, i));
    
    if ~exist(hrfname, 'file')
        disp('You must have already written HRF images in the current dir.  See fmri_model.robustfit.');
    end
    
    hrfs{i} = fmri_data(hrfname, mask);
    
end

% voxels
v = size(hrfs{1}.X, 2);

% test build to set up sizes
for ii = 1:nconds
    voxelhrfs{ii} = hrfs{ii}.X(:, i);
end
obj = build_single_trial(obj, voxelhrfs);
nb = length(obj.xX.iH);
betas = zeros(nb, v, 'single');


fprintf('Fitting single-trial for %3.0f voxels: 000000', v);
t1 = tic;

for i = 1:v
    
    if mod(i, 100) == 0 || i == v
        fprintf('\b\b\b\b\b\b%06d', i);
    end
    
    for ii = 1:nconds
        voxelhrfs{ii} = hrfs{ii}.X(:, i);
    end
    
    % Build model
    % ---------------------------------------------------------------------
    obj = build_single_trial(obj, voxelhrfs);
    
    
    % Fit OLS model (trial regressors to data)
    % ---------------------------------------------------------------------
    
    b = pinv(obj.xX.X) * fmri_data_obj.X(:, i);
    betas(:, v) = b(1:nb);
end

t2 = toc(t1);
fprintf(' Done in %3.0f s\n', t2);

% Parse betas into conditions and write images
% ---------------------------------------------------------------------
c = obj.xX.cond_assignments;
nconds = size(c, 2); %***need to add pms and use above
vI = mask.volInfo;

for i = 1:nconds
    imgname = sprintf('single_trial_amp_cond%3.0f_%s.img', i, obj.xX.name{i});
    
    write_image(b(c(:, i), :)', imgname, vI, badvox);
end



end % function