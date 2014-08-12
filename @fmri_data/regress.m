function [out, statimg] = regress(dat, varargin)
% [out, statimg] = regress(dat, [p-val threshold], [thresh_type], ['nodisplay'])
%
% regression method for fmri_data object
% regress dat.Y on dat.dat at each voxel, and return voxel-wise statistic
% images.
%
% Each column of Y is a predictor in a multiple regression,
% and the intercept is the last column.
% Thus, the output images correspond to the columns of Y, with the last
% column being the intercept.
% The reason the "Y" field is used as the "X" design matrix is that Y in
% the fmri_data object is reserved for outcome data/behavior.
% In the standard brain-mapping approach, Y is the predictor and each
% voxel's data is the outcome in a series of regressions.
%
% This function can also create a map of brain regions that predict the Y
% vector using the 'brainony' option.  This is essentially a univariate
% version of the 'predict' command.  There is also a robust option for this
% ('brainony', 'robust')
%
% Inputs:
% dat should be an fmri_data object with Y field defined.
% optional inputs in [  ] above are:
% p-value threshold
% string indicating threshold type (see help statistic_image.threshold for options)
%
% Optional Inputs:
% 'nointercept'     estimates model without intercept
% 'nodisplay'       does not display resulting model
% 'brainony'        univariate approach to predict obj.Y from brain data
%
% Outputs:
% out is a structure of information about the regression
% statimg is a statistic_image object that can be thresholded and
% plotted/imaged.  statimg.dat contains T-values, .p contains p-values for
% each regressor.
%
% Examples:
% % Run regression with liberal threshold
% [out, statimg] = regress(data_comb, .05, 'unc');
%
% % Re-threshold at different values
% statimg = threshold(statimg, .05, 'fdr');
% statimg = threshold(statimg, .001, 'unc');
%
% % Re-display results of thresholding
% orthviews(statimg);
%
% %Run a regression predicting behavior from brain at liberal threshold
% [out, statimg]  = regress(data_comb, .05, 'unc', 'brainony')
%
% %Run a robust regression predicting behavior from brain at liberal
% threshold.  WARNING! Very slow!
% [out, statimg]  = regress(data_comb, .05, 'unc', 'brainony','robust')

%Notes:
% c Tor Wager, Dec 2010
% Edited by Luke Chang, 9/27/2012 to add optional input to reverse X & Y (i.e., create a map of voxels that predict the behavioral variable)
% Edited by Luke Chang, 9/28/2012 to add optional input to run robust regression for brainony
% Edited by Luke Chang, 10/24/2012 to save residuals (i.e., out.r), which is helpful for denoising an image
% Edited by Luke Chang, 3/26/2013 to add optional input to not add an intercept - allows for more flexible modeling options


% default options for thresholding
inputargs = {[.001], 'uncorrected'};

dodisplay = 1;

dobrainony = 0; %added by lc
dorobust = 0; %added by lc
dointercept = 1; %added by lc

if ~isempty(varargin)
    for i = 1:length(varargin)
        inputargs{i} = varargin{i};
    end
end

if any(strcmp(inputargs, 'nodisplay'))
    i = find(strcmp(inputargs, 'nodisplay'));
    dodisplay = 0;
    inputargs(i) = [];
end

if any(strcmp(inputargs, 'brainony'))
    i = find(strcmp(inputargs, 'brainony'));
    dobrainony = 1;
    inputargs(i) = [];
end

if any(strcmp(inputargs, 'robust'))
    i = find(strcmp(inputargs, 'robust'));
    dorobust = 1;
    inputargs(i) = [];
end

if any(strcmp(inputargs, 'nointercept'))
    i = find(strcmp(inputargs, 'nointercept'));
    dointercept = 0;
    inputargs(i) = [];
end

n = size(dat.dat, 2);
if n ~= size(dat.Y, 1), error('dat.dat must have same number of columns as dat.Y has rows.'); end

if ~dobrainony %default is to regress Y on brain
    
    if dointercept  %can fit model with no intercept: 3/26/13 LC
        X = dat.Y; X(:, end+1) = 1;
    else
        X = dat.Y;
    end
    
    % display info about regression
    fprintf('regression > Y: %3.0f voxels. X: %3.0f obs, %3.0f regressors, intercept is last.\n', size(dat.dat, 1), n, size(X, 2));
    
    % model fit
    tic
    [n, k] = size(X);
    b = pinv(X) * dat.dat';
    
    % error
    r = dat.dat' - X*b;
    sigma = std(r);
    
    stderr = ( diag(inv(X' * X)) .^ .5 ) * sigma;  % params x voxels matrix of std. errors
    
    % inference
    
    t = b ./ stderr;
    
    dfe = n - k;
    p = 2 * (1 - tcdf(abs(t), dfe));
    
    sigma(sigma == 0) = Inf;
    t(isnan(t)) = 0;
    p(isnan(p)) = 0;
    
    toc
else %regress brain on Y
    if ~dorobust %Run standard OLS regression
        % display info about regression
        fprintf('Predicting dat.Y from Brain Activity');
        
        tic
        warning off
        for  i=1:size(dat.dat,1) %slow, need to figure out how to vectorize this
            
            if dointercept  %can fit model with no intercept: 3/26/13 LC
                X = dat.Y; X(:, end+1) = 1;
            else
                X = dat.Y;
            end
            
            % model fit
            b(:,i) = pinv(X)*dat.Y;
            
            %error
            r(:,i)=dat.Y-X*b(:,i);
            sigma(i)=std(r(:,i));
            stderr(:,i) = (diag(inv(X'*X)).^ .5) * sigma(i); %warning about matrix being singular not sure what the problem is.
        end
        warning on
        [n, k] = size(X);
        % inference
        dfe = n - k;
        
        t = b ./ stderr;
        p = 2 * (1 - tcdf(abs(t), dfe));
        
        sigma(sigma == 0) = Inf;
        t(isnan(t)) = 0;
        p(isnan(p)) = 0;
        
        toc
        
    else  %Run robust regression %Warning!  Very Slow using the loop
        
        fprintf('Predicting dat.Y from Brain Activity Using Robust Regression');
        tic
        warning off
        for  i=1:size(dat.dat,1) %slow, need to figure out how to vectorize this
            X=dat.dat(i,:)'; X(:, end+1) = 1;
            Y=dat.Y;
            
            % model fit
            [bb,stats]=robustfit(X, Y, 'bisquare', [], 'off');
            b(:,i)=bb;
            t(:,i)=stats.t;
            p(:,i)=stats.p;
            stderr(:,i)=stats.se;
            sigma(:,i)=stats.robust_s; %robust estimate of sigma. LC not sure this is the best choice can switch to 'OLS_s','MAD_s', or 's'
            
        end
        warning on
        
        [n, k] = size(X);
        dfe = stats.dfe;
        
        sigma(sigma == 0) = Inf;
        t(isnan(t)) = 0;
        p(isnan(p)) = 0;
        
        toc
        
    end
end

% save results
statimg = statistic_image;
statimg.type = 'T';
statimg.p = p';
statimg.ste = stderr';
statimg.N = n;
statimg.dat = t';
statimg.dat_descrip = sprintf('t-values from OLS regression, intercept is last');
statimg.volInfo = dat.volInfo;
statimg.removed_voxels = dat.removed_voxels;
statimg.removed_images = false;  % this image does not have the same dims as the original dataset

statimg = threshold(statimg, inputargs{:});

out.X = X;
out.b = b;
out.sigma = sigma;
out.stderr = stderr;
out.dfe = dfe;
out.t = t;
out.p = p;
out.r = r; %LC : 10/24/12 : save residual - helpful for denoising


% view results
if k < 10 & dodisplay
    
    orthviews(statimg);
    spm_orthviews_name_axis('Intercept', size(statimg.dat, 2))
    
end

end % function