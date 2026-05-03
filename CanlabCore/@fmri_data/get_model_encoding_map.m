function [beta,p_values] = get_model_encoding_map(roi, latent_timeseries_pathway)
% get_model_encoding_map Voxelwise encoding map relating ROI activity to MPathI latent time series.
%
% :Usage:
% ::
%
%     [beta, p_values] = get_model_encoding_map(roi, latent_timeseries_pathway)
%
% Computes a model encoding map by performing linear regression to
% understand how the latent time-series data from MPathI
% (model_brain_pathway.m) relate to voxel activity data (ROI). For
% each voxel and each pathway, beta coefficients and p-values are
% returned, revealing how fluctuations in the latent time series impact
% neural activity patterns and forming the model encoding map.
%
% Each beta quantifies the change in activity at a given voxel
% associated with a one-unit change in the corresponding latent time
% series, i.e., the strength and direction of the linear relationship
% between voxel activity and the latent pathway.
%
% :Inputs:
%
%   **roi:**
%        Actual activity time-series data, shape [TR x voxels].
%
%   **latent_timeseries_pathway:**
%        Latent time-series data from MPathI (XS or YS), shape
%        [TR x pathways] (e.g., two columns representing on- and
%        off-target conditions).
%
% :Outputs:
%
%   **beta:**
%        Beta coefficients describing the linear relationship between
%        each voxel's activity and each latent time series pathway,
%        shape [voxels x pathways].
%
%   **p_values:**
%        p-values for the significance of each beta coefficient,
%        indicating the strength of each relationship.
%
% :Examples:
% ::
%
%     stats = model_brain_pathway(preproc_dat, ...
%         source_one, source_two, target_one, target_two, ...
%         'Indices', wh_run * (1:size(wh_run, 2))');
%
%     score{1} = [stats.latent_timeseries_source(:,1), stats.latent_timeseries_source(:,3)]; % source1
%     score{2} = [stats.latent_timeseries_target(:,1), stats.latent_timeseries_target(:,2)]; % target1
%     score{3} = [stats.latent_timeseries_source(:,4), stats.latent_timeseries_source(:,2)]; % source2
%     score{4} = [stats.latent_timeseries_target(:,4), stats.latent_timeseries_target(:,3)]; % target2
%
%     roi{1} = stats.source_one_obj.dat';
%     roi{2} = stats.target_one_obj.dat';
%     roi{3} = stats.source_two_obj.dat';
%     roi{4} = stats.target_two_obj.dat';
%
%     for k = 1:4   % S1, T1, S2, T2
%         [b_temp, p_temp] = get_model_encoding_map(roi{k}, score{k});
%     end
%
% :See also:
%   - model_brain_pathway, model_mpathi
%
% ..
%    Developed by Tor and Byeol, 2022-24.
% ..
% Get the dimensions of the input data
[tr_n, voxels_n] = size(roi);                   % [TR x voxels]
pathway_n = size(latent_timeseries_pathway,2);  % [TR x pathways]

I  = ones(tr_n, 1); % intercept
beta = zeros(voxels_n, pathway_n);      
p_values = zeros(voxels_n, pathway_n);  

% Perform linear regression for each pathway and voxel
for j = 1:pathway_n
    b = zeros(2, voxels_n);
    for i = 1:voxels_n
        % Perform linear regression
        X = [roi(:, i) I]; % one voxel + intercept
        Y = latent_timeseries_pathway(:, j);

        [b(:,i), ~,~,~, stats] = regress(Y, X);
        
        % Calculate p-value for the beta coefficient
        p_values(i,j) = stats(3);
    end
    beta(:,j) = b(1,:)'; 
end

end