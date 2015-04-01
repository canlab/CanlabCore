function out = regress(dat, varargin)
% out = regress(dat, varargin)
%
% regression method for fmri_data object
%
% Regress dat.X on dat.dat at each voxel, and return voxel-wise statistic
% images. Each column of dat.X is a predictor in a multiple regression,
% and the intercept is the last column. Intercept will automatically be
% added if not detected unless 'nointercept' is specified.  
%
% This function can also create a map of brain regions that predict the dat.Y
% vector using the 'brainony' option.  This is essentially a univariate
% version of the 'predict' command.  Warning: this is very slow as it loops
% through all voxels.
%
% Regression is OLS by default, but can be robust using 'robust' flag.
% Warning: Running robust option is considerably slower than OLS.
%
% Creates thresholded plot by default
% ---------------------------------------------------------------------
% Inputs:
% ---------------------------------------------------------------------
% dat:                  should be an fmri_data object with X field defined.
%                       dat.X can be a design_matrix() object.
%
% ---------------------------------------------------------------------
% Optional Inputs:
% ---------------------------------------------------------------------
% [threshold, 'unc']:   p-value threshold string indicating threshold type
%                       (see help statistic_image.threshold for options)
% nointercept:          Do not add intercept to model
% nodisplay:            Do not plot thresholded results using orthviews
% brainony:             univariate approach to predict obj.Y from brain data
% residual:             Output residual as fmri_data() object
%
% ---------------------------------------------------------------------
% Outputs:
% ---------------------------------------------------------------------
% out:                  A structure containing stats_img and fmri_data
%                       objects.
% out.b:                stats_img object of beta values estimated from
%                       regression
% out.t:                stats_img object of t-values with input threshold
% out.df:               fmri_data object of degrees of freedom
% out.sigma:            fmri_data object of variance of residual
% out.residual:         (optional) fmri_data object of residual data after
%                       model has been regressed out.
%
% ---------------------------------------------------------------------
% Examples:
% ---------------------------------------------------------------------
% % Run regression with liberal threshold
% out = regress(dat, .05, 'unc');
%
% % Run regression with conservative threshold and save residual
% out = regress(dat, .001, 'unc', 'residual);
%
% % Run robust regression with fdr threshold
% out = regress(dat, .05, 'fdr','robust');
%
% %Run a regression predicting behavior from brain at liberal threshold
% out  = regress(data_comb, .05, 'unc', 'brainony')
%
% % Re-threshold at different values
% out.t = threshold(out.t, .05, 'fdr');
% out.t = threshold(out.t, .001, 'unc');
%
% % Re-display results of thresholding
% orthviews(out.t);
%
% % Write out beta image to current directory
% out.b.fullpath = fullfile(pwd,'beta.nii');
% write(out)
%
% ---------------------------------------------------------------------
% Copyright (c) 2015 Tor Wager & Luke Chang
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Notes:
% c Tor Wager, Dec 2010
% Edited by Luke Chang, 9/27/2012 to add optional input to reverse X & Y (i.e., create a map of voxels that predict the behavioral variable)
% Edited by Luke Chang, 9/28/2012 to add optional input to run robust regression for brainony
% Edited by Luke Chang, 10/24/2012 to save residuals (i.e., out.r), which is helpful for denoising an image
% Edited by Luke Chang, 3/26/2013 to add optional input to not add an intercept - allows for more flexible modeling options
% Code completely refactored by Luke Chang 2/24/25

% ---------------------------------------------------------------------
% Defaults
% ---------------------------------------------------------------------
inputargs = {[.001], 'uncorrected'}; % default options for thresholding
do_display = 1;
do_x_on_brain = 1; %else do brain on Y
do_robust = 0;
do_intercept = 1;
do_resid = 0;

% ---------------------------------------------------------------------
% Parse Inputs
% ---------------------------------------------------------------------
for varg = 1:length(varargin)
    if ischar(varargin{varg})
        % reserved keywords
        if strcmpi('unc',varargin{varg})
            inputargs = {varargin{varg-1}, 'uncorrected'};
            varargin{varg} = {}; varargin{varg - 1} = {};
        end
        if strcmpi('fdr',varargin{varg})
            inputargs = {varargin{varg-1}, 'fdr'};
            varargin{varg} = {}; varargin{varg - 1} = {};
        end
        if strcmpi('nodisplay',varargin{varg})
            do_display = 0;
            varargin{varg} = {};
        end
        if strcmpi('robust',varargin{varg})
            do_robust = 1;
            varargin{varg} = {};
        end
        if strcmpi('brainony',varargin{varg})
            do_x_on_brain = 0;
            varargin{varg} = {};
        end
        if strcmpi('nointercept',varargin{varg})
            do_intercept = 0;
            varargin{varg} = {};
        end
        if strcmpi('residual',varargin{varg})
            do_resid = 1;
            varargin{varg} = {};
        end
    end
end

% ---------------------------------------------------------------------
% Check Data
% ---------------------------------------------------------------------

% Check if fmri_data or image_vector
if ~isa(dat,'fmri_data')
    error('dat input must be fmri_data object')
end

% Check of dat.X exists and is correct format
if do_x_on_brain
    if isempty(dat.X)
        error('Make sure you include a design matrix in dat.X')
    end
    if size(dat.dat, 2) ~= size(dat.X, 1)
        error('dat.dat must have same number of columns as dat.X has rows.')
    end
    if isa(dat.X,'design_matrix')
       dat.X = dat.X.dat; 
    end
else % Check if dat.Y exists and is in correct format if running brainony
    if isempty(dat.Y)
        error('Make sure you include a vector in dat.Y.')
    end
    if size(dat.dat, 2) ~= size(dat.Y, 1)
        error('dat.dat must have same number of columns as dat.Y has rows.')
    end
end

% Check if Rank Deficient
if rank(dat.X) < size(dat.X,2)
    sprintf('Warning:  dat.X is rank deficient.')
end

% Check if Intercept is in model or specified for x_on_brain default
if do_intercept && do_x_on_brain
    wh_int = intercept(dat.X,'which');
    
    if isempty(wh_int)
        sprintf('No intercept detected, adding intercept to last column of design matrix')
        X = intercept(dat.X,'add');
    else
        sprintf('Intercept detected in column %1.0f of dat.x',wh_int)
        X = dat.X;
    end
    
else
    X = dat.X;
end

% ---------------------------------------------------------------------
% Run Regression
% ---------------------------------------------------------------------

tic
warning off

if do_x_on_brain %default is to regress dat.X on dat.dat (x on brain)
    % display info about regression
    fprintf('regression > Y: %3.0f voxels. X: %3.0f obs, %3.0f regressors, intercept is last.\n', size(dat.dat, 1), size(X, 2));
    fprintf('\nPredicting Brain Activity from dat.X');
    
    if do_robust %need to loop through voxels - Slow!
        fprintf('\nRunning in Robust Mode');
        
        for i = 1:size(dat.dat,1)
            
            [n, k] = size(X);
            [bb,stats] = robustfit(X, dat.dat(i,:)', 'bisquare', [], 'off');
            
            b(:,i)=bb; %Betas
            t(:,i)=stats.t; %t-values
            p(:,i)=stats.p; %p-values
            dfe(:,i)=stats.dfe; %degrees of freedom
            stderr(:,i)=stats.se; %Standard Error
            sigma(:,i)=stats.robust_s; %robust estimate of sigma. LC not sure this is the best choice can switch to 'OLS_s','MAD_s', or 's'
        end
        r = dat.dat' - X*b; %residual
        
    else %OLS - vectorized - Fast!
        fprintf('\nRunning in OLS Mode');
        
        % Estimate Betas in vector
        [n, k] = size(X);
        b = pinv(X) * dat.dat';
        
        % Error
        r = dat.dat' - X*b;
        sigma = std(r);
        stderr = ( diag(inv(X' * X)) .^ .5 ) * sigma;  % params x voxels matrix of std. errors
        
        % Inference
        [t,dfe,p,sigma] = param_t_test(X,b,stderr,sigma);
        
    end
else % Regress brain on Y - loops through voxels, slow!
    
    % display info about regression
    fprintf('regression > X: %3.0f voxels. Y: %3.0f obs, %3.0f regressors, intercept is last.\n', size(dat.dat, 1), size(dat.Y, 2));
    fprintf('\nPredicting dat.Y from Brain Activity');
    
    if do_robust %need to loop through voxels - Slow!
        fprintf('\nRunning in Robust Mode');
        for i = 1:size(dat.dat,1)
            
            % Create X from brain Data
            if do_intercept
                X = intercept(dat.dat(i,:)','add');
            else
                X = dat.dat(i,:)';
            end
            [n, k] = size(X);
            [bb,stats] = robustfit(X, dat.Y, 'bisquare', [], 'off');
            
            b(:,i)=bb; %Betas
            t(:,i)=stats.t; %t-values
            p(:,i)=stats.p; %p-values
            dfe(:,i)=stats.dfe; %degrees of freedom
            stderr(:,i)=stats.se; %Standard Error
            sigma(:,i)=stats.robust_s; %robust estimate of sigma. LC not sure this is the best choice can switch to 'OLS_s','MAD_s', or 's'
            r(:,i) = dat.Y - X * b(:,i); %residual
        end
        
    else %OLS
        fprintf('\nRunning in OLS Mode');
        
        for i = 1:size(dat.dat,1)
            
            % Create X from brain Data
            if do_intercept
                X = intercept(dat.dat(i,:)','add');
            else
                X = dat.dat(i,:)';
            end
            
            % Estimate Betas in vector
            b(:,i) = pinv(X) * dat.Y;
            
            % Error
            r(:,i) = dat.Y - X * b(:,i);
            sigma(i) = std(r(:,i));
            stderr(:,i) = ( diag(inv(X' * X)) .^ .5 ) * sigma(i);  % params x voxels matrix of std. errors
        end
        % Inference
        [t,dfe,p,sigma] = param_t_test(X,b,stderr,sigma);
    end
end
stop = toc;
fprintf('\nModel run in %d minutes and %.2f seconds\n',floor(stop/60),rem(stop,60));

% ---------------------------------------------------------------------
% Create Output
% ---------------------------------------------------------------------

out = struct;

% Betas
out.b = statistic_image;
out.b.type = 'Beta';
out.b.p = p';
out.b.ste = stderr';
out.b.N = n;
out.b.dat = b';
out.b.dat_descrip = sprintf('Beta Values from regression, intercept is last');
out.b.volInfo = dat.volInfo;
out.b.removed_voxels = dat.removed_voxels;
out.b.removed_images = false;  % this image does not have the same dims as the original dataset
out.b = threshold(out.b, inputargs{:}); % Threshold image

% T stats
out.t = statistic_image;
out.t.type = 'T';
out.t.p = p';
out.t.ste = stderr';
out.t.N = n;
out.t.dat = t';
out.t.dat_descrip = sprintf('t-values from regression, intercept is last');
out.t.volInfo = dat.volInfo;
out.t.removed_voxels = dat.removed_voxels;
out.t.removed_images = false;  % this image does not have the same dims as the original dataset
out.t = threshold(out.t, inputargs{:}); %Threshold image

% DF as fmri_data
out.df = dat;
out.df.dat = dfe';
out.df.dat_descrip = sprintf('Degrees of Freedom');

% Sigma as fmri_data
out.sigma = dat;
out.sigma.dat = sigma';
out.sigma.dat_descrip = sprintf('Sigma from Regression');

% Residual as fmri_data
if do_resid
    out.resid = dat;
    out.resid.dat = r';
    out.resid.dat_descrip = sprintf('Residual from Regression');
end

% ---------------------------------------------------------------------
% Plot Results
% ---------------------------------------------------------------------
if k < 10 && do_display
    
    orthviews(out.t);
    
end

warning on

% ---------------------------------------------------------------------
% Subfunctions
% ---------------------------------------------------------------------

    function [t,dfe,p,sigma] = param_t_test(X,b,stderr,sigma)
        % test whether parameter is significantly different from zero
        %
        % Inputs:
        % X:        design matrix
        % b:        beta values
        % stderr:   standard error of beta estimate
        % sigma:    standard deviation of residual
        %
        % Returns:
        % t:        t-value
        % dfe:      degrees of freedom
        % p:        p-value
        % sigma:    standard deviation of residual
        
        [n, k] = size(X);
        t = b ./ stderr;
        dfe = n - k;
        p = 2 * (1 - tcdf(abs(t), dfe));
        
        sigma(sigma == 0) = Inf;
        t(isnan(t)) = 0;
        p(isnan(p)) = 0;
        dfe = repmat(dfe,1,size(t,2));
    end

end


