function t_obj = parcel_stats2statistic_image(parcel_obj, tscores, pvalues, dfe, sig)
% Utility that transforms parcel-wise t-statistics and p-values into voxelwise image space
%
% t_obj = parcel_stats2statistic_image(atlas_obj, input_struct)
%
%parcel_obj = PARCELS.Shen.parcel_obj;
% t_obj = parcel_stats2statistic_image(parcel_obj, input_struct );

% Create placeholder statistic image

parcel_obj = replace_empty(parcel_obj);   
k = size(tscores, 2);
placeholder_vec = ones(parcel_obj.volInfo.n_inmask, k);
all_parcel_idx = double(parcel_obj.dat);
u = unique(all_parcel_idx); u(u == 0) = [];

% initialize variables with correct size

t_obj = statistic_image('dat', single(0 .* placeholder_vec), ...
    'p', placeholder_vec, ...
    'sig', logical(placeholder_vec), ...
    'type', 'T', ...
    'dfe', placeholder_vec, ...
    'volInfo', parcel_obj.volInfo);

% number of parcels 

n = size(tscores, 1); % num parcels

if n ~= length(u)
        error('Parcel indices in parcel_obj and extracted parcel-wise t-values do not match. Check code.')
end



t_obj.dat = zeros(size(placeholder_vec, 1), k); % voxels x conditions


    for j = 1:length(u)
        % For each parcel, fill in statistic_image object
        % ----------------------------------------------------
        parcelidx = u(j);
        
        wh_vox = all_parcel_idx == parcelidx;
        
        % map parcels to voxels
        t_obj.dat(wh_vox, :) = repmat(tscores(j, :), sum(wh_vox), 1);
        t_obj.p(wh_vox, :) = repmat(pvalues(j, :), sum(wh_vox), 1);
        t_obj.dfe(wh_vox, :) = repmat(dfe(j, 1), sum(wh_vox), k);       % replicate dfe
        t_obj.sig(wh_vox, :) = repmat(sig(j, :), sum(wh_vox), 1);
                
    end
       
t_obj = enforce_variable_types(t_obj);  % space-saving: 5/24/17

% input_struct.t_statistic_obj = t_obj;

end % function
