
function [spatial_maps, timecourses, tmaps, GM_dat] = dual_regression_fsl(GM, T, varargin)
% DUAL_REGRESSION_FSL Perform FSL-style dual regression on fMRI data
%
% :Usage:
% ::
%    [spatial_maps, timecourses, tmaps] = dual_regression_fsl(GM, T, 'zscore_data', true, 'n_iter', 10, 'add_intercept', false);
%
% :Inputs:
%   **GM** : [voxels x components] Group ICA spatial maps, or an fmri_data object
%   **T**  : [voxels x timepoints] Subject fMRI data, or an fmri_data object
%
% :Optional Inputs:
%   **'zscore_data'**   : [logical] Whether to z-score timecourses and spatial maps. Default = true.
%   **'zscore_maps'**   : [logical] Whether to z-score timecourses and spatial maps. Default = true.
%   **'n_iter'**        : [integer >= 1] Number of dual regression iterations. Default = 1.
%   **'add_intercept'** : [logical] Whether to add an intercept (column of ones) to GM. Default = true.
%   **'verbose'**       : [logical] Verbose output. Default = true.
%   **'doplot'**        : [logical] Plot maps and timecourses. Default = true.
%
% :Outputs:
%   **spatial_maps** : Subject-specific spatial maps [voxels x components] or fmri_data
%   **timecourses**  : Subject-specific component timecourses [components x timepoints]
%   **tmaps**        : z-scored spatial maps; same format as spatial_maps
%

% Detect fmri_data objects
is_fmri_GM = isa(GM, 'fmri_data');
is_fmri_T  = isa(T, 'fmri_data');

if is_fmri_GM
    GM_dat = double(GM.dat);
else
    GM_dat = double(GM);
end

if is_fmri_T
    T_dat = double(T.dat);
else
    T_dat = double(T);
end

% Check matrix dimensions
if size(GM_dat, 1) ~= size(T_dat, 1)
    error('GM and T must have the same number of voxels (rows).');
end

% Parse inputs
p = inputParser;
addParameter(p, 'zscore_data', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'zscore_maps', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'n_iter', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'add_intercept', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'verbose', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'doplot', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
zscore_data = p.Results.zscore_data;
zscore_maps = p.Results.zscore_data;
n_iter = p.Results.n_iter;
add_intercept = p.Results.add_intercept;
verbose = p.Results.verbose;
doplot = p.Results.doplot;

% Z-score spatial maps
if zscore_maps
    if verbose, fprintf('Z-scoring group spatial maps...\n'); end
    GM = zscore(GM, 0, 1);
end

% Mean-centering or intercept
    if verbose, fprintf('Adding intercept to group maps...\n'); end
    GM_orig_dat = GM_dat; % after any scaling/etc
    GM_dat = [GM_dat, ones(size(GM_dat, 1), 1)];

if p.Results.zscore_data  % GM_orig_dat for plotting only, for comparability 
    GM_orig_dat = zscore(GM_orig_dat);
end

% ----------------------------------------------------------------- 
% Iterations
% -----------------------------------------------------------------

for iter = 1:n_iter
    if verbose
        fprintf('Dual regression iteration %d/%d\n', iter, n_iter);
    end

    % Step 1: Spatial regression
    if verbose, fprintf('   Spatial regression of data on spatial maps\n'); end
    timecourses = pinv(GM_dat) * T_dat;  % maps x time (data series)


    % Prepare time courses as design matrix for next step
    if zscore_data
        if verbose, fprintf('   Z-scoring each timecourse...\n'); end
        timecourses = zscore(timecourses, 0, 2);

    elseif add_intercept
        timecourses = [timecourses(1:end-1, :); ones(1, size(timecourses,1))];  

    else
        warning('You should either z-score data or add intercept!')
    end


    % Step 2: Temporal regression
    % 
    % This is the same as spatial_maps = (pinv(timecourses') * T_dat')'
    % So a regression of each voxel's data onto each timecourse, but
    % transposed

    if verbose, fprintf('   Temporal regression of data onto time courses...\n'); end
    spatial_maps = T_dat * timecourses' / (timecourses * timecourses');

    % if add_intercept
    %     spatial_maps = spatial_maps(1:end - 1, :); % check
    % end

    % spatial_maps2 = (pinv(timecourses') * T_dat')';
    % d = spatial_maps - spatial_maps2; max(abs(d(:)))

    % Update GM for next iteration (excluding intercept column if used)

    if zscore_maps
        if verbose, fprintf('   Z-scoring spatial maps...\n'); end
        spatial_maps = zscore(spatial_maps, 0, 1);

    elseif add_intercept
        GM_dat = [spatial_maps(:, 1:end-1), ones(size(spatial_maps,1),1)];

    else
        GM_dat = spatial_maps;
        warning('You should either z-score maps or add intercept!')
    end


end % iterations

% final regression to get final timecourses
timecourses = pinv(GM_dat) * T_dat;  % maps x time (data series)

if add_intercept
    spatial_maps = spatial_maps(:, 1:end-1);
    timecourses = [timecourses(1:end-1, :)];
end

% **** t-maps
tmaps = spatial_maps;

% Wrap in fmri_data if needed
if is_fmri_T
    spatial_maps = fmri_data(spatial_maps, T.volInfo, T.mask);
    tmaps = fmri_data(tmaps, T.volInfo, T.mask);
end

% Optional plotting
if doplot

    if is_fmri_T
        spatial_maps_dat = spatial_maps.dat;
    else
        spatial_maps_dat = spatial_maps;
    end

    create_figure('Spatial Maps Comparison', 1, 3);
    subplot(1,3,1); imagesc(GM_orig_dat); title('Original Group Maps'); axis tight
    xlabel('Components'); ylabel('Voxels'); colorbar;

    subplot(1,3,2); imagesc(spatial_maps_dat); title('Reconstructed Spatial Maps');  axis tight
    xlabel('Components'); ylabel('Voxels'); colorbar;
    subplot(1,3,3); scatter(GM_orig_dat(:), spatial_maps_dat(:), 10, '.');
    r = corr(GM_orig_dat(:), spatial_maps_dat(:));
    title(sprintf('Spatial Regression, r = %3.2f'), r)
    xlabel('Original Maps'); ylabel('Reconstructed');

    projected_tc = zscore(GM_orig_dat' * T_dat);
    create_figure('Time Course Comparison', 1, 3);
    subplot(1,3,1); plot_matrix_cols(projected_tc'); title('Projected Timecourses (GM^T * T)');
    subplot(1,3,2); plot_matrix_cols(timecourses'); title('Estimated Timecourses');
    subplot(1,3,3); scatter(projected_tc(:), timecourses(:), 10, '.');
    r = corr(projected_tc(:), timecourses(:));
    title(sprintf('Temporal Regression, r = %3.2f'), r)
    xlabel('Projected'); ylabel('Estimated');

    % d = GM_orig_dat - spatial_maps_dat;
    % max(abs(d(:)))
    % figure; imagesc(d); colorbar

end

end % function

