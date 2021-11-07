function [est_outliers_uncorr, est_outliers_corr, outlier_table, outlier_regressor_matrix_uncorr, outlier_regressor_matrix_corr] = outliers(dat, varargin)
% Find outliers based on a combination of rmssd (DVARS), robust spatial
% variance (MAD), Mahalanobis distance based on covariance and correlation
% matrices, and images with >25% missing values
%
% :Usage:
% ::
%
% [est_outliers_uncorr, est_outliers_corr, outlier_table, outlier_regressor_matrix_uncorr, outlier_regressor_matrix_corr] = outliers(dat, ['noverbose'])
%
% Images usually change slowly over time, and sudden changes in intensity can also often be a sign of bad things
% -- head movement artifact or gradient misfires, interacting with the magnetic field to create distortion
% across the brain.
%
% RMSSD (DVARS) tracks large changes across successive images, regardless of what the sign of the changes is or where they are.
% In addition, images with unusually high spatial standard deviation across voxels may be outliers with image
% intensity distortions in some areas of the image but not others (e.g., bottom half of brain vs. top half,
% or odd vs. even slices).
% This matrix can be added to your design matrix as a set of nuisance covariates of no interest.
%
% *Optional Inputs:
%
%   **'madlim'**
%   Limit for median absolute deviation (MAD) for rmssd and spatial abs. deviation
%
% :Examples:
% ::
% obj = load_image_set('kragel18_alldata');
% obj2 = rescale(obj, 'l2norm_images');     % normalize heterogeneous datasets
% [est_outliers_uncorr, est_outliers_corr, outlier_table, outlier_regressor_matrix_uncorr, outlier_regressor_matrix_corr] = outliers(obj);
%

% Programmers Notes:
% ..
%    Created 11/6/2021 by Tor Wager, from a combination of other code
%    (default: based on estimated outliers at 3 standard deviations.)
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------
madlim = 3;

doverbose = true;
verbosestr = 'doverbose';

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'madlim' 'doverbose'};

keyword_inputs = {'noverbose'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            case keyword_inputs
                % Skip, deal with these below
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% 2nd pass: Keyword inputs. These supersede earlier inputs
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
      
            case {'noverbose'}
                doverbose = false;
                
        end
    end
end

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

if doverbose
    fprintf('rmssd: ');
else
    verbosestr = 'noverbose';
end

% Get RMSSD (DVARS)
sdiffs = diff(dat.dat')';
sdiffs = [mean(sdiffs, 2) sdiffs]; % keep in image order

rmssd = ( mean(sdiffs .^ 2) ) .^ .5;    % rmssd - root mean square successive diffs
rmssd(1) = median(rmssd);               % avoid first time point being very different and influencing distribution and plots.

rmssd_outliers = rmssd' > mean(rmssd) + madlim * mad(rmssd);

if doverbose
    fprintf('%3.0f images ', sum(rmssd_outliers));
end

if doverbose
    fprintf('spatial mad: ');
end
%mymean = mean(dat.dat); % global mean
spatialmad = mad(dat.dat);  % variation in successive diffs - like rmssd but without abs mean diff. artifacts in some parts of image...
spatialmad(1) = median(spatialmad);
spatialmad_outliers = spatialmad' > mean(spatialmad) + madlim * mad(spatialmad);

if doverbose
    fprintf('%3.0f images ', sum(spatialmad_outliers));
end


if doverbose
    fprintf('Missing values: ');
end

desc = descriptives(dat, 'noverbose');
missingvals = desc.images_missing_over_25percent';

if doverbose
    fprintf('%3.0f images ', sum(missingvals));
    fprintf('\n\n');
end

if doverbose
    fprintf('Mahalanobis (cov and corr, q<0.05 corrected):\n');
end
% Get mahalanobis outliers
[mahalcov, ~, ~, mahal_cov_outlier_uncorr, mahal_cov_outlier_corr] = mahal(dat, 'noplot', verbosestr);
[mahalcorr, ~, ~, mahal_corr_outlier_uncorr, mahal_corr_outlier_corr] = mahal(dat, 'corr', 'noplot', verbosestr);

if doverbose
    fprintf('%3.0f images ', sum(mahal_cov_outlier_corr | mahal_corr_outlier_corr));
end


est_outliers_uncorr = rmssd_outliers | spatialmad_outliers | mahal_cov_outlier_uncorr | mahal_corr_outlier_uncorr | missingvals;
est_outliers_corr = rmssd_outliers | spatialmad_outliers | mahal_cov_outlier_corr | mahal_corr_outlier_corr | missingvals;

% Make Table

outlier_table = table(est_outliers_uncorr, est_outliers_corr, rmssd', spatialmad', mahalcov, mahalcorr, rmssd_outliers, spatialmad_outliers, mahal_cov_outlier_uncorr, mahal_cov_outlier_corr, mahal_corr_outlier_uncorr, mahal_corr_outlier_corr, ...
    'VariableNames', {'outliers_uncorr' 'outliers_corr' 'rmssd_dvars' 'spatialmad' 'mahal_cov' 'mahal_corr' 'rmssd_dvars_indic', 'spatialmad_indic', 'mahal_cov_uncor_indic', 'mahal_cov_cor_indic', 'mahal_corr_uncor_indic', 'mahal_corr_cor_indic'});
   
% Make indicator matrices

outlier_regressor_matrix_uncorr = intercept_model(size(dat.dat, 2), find(est_outliers_uncorr));
outlier_regressor_matrix_uncorr = outlier_regressor_matrix_uncorr(:, 2:end);  % remove initial intercept column

outlier_regressor_matrix_corr = intercept_model(size(dat.dat, 2), find(est_outliers_corr));
outlier_regressor_matrix_corr = outlier_regressor_matrix_corr(:, 2:end);  % remove initial intercept column

end % main function


