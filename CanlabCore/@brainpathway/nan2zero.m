function bs = nan2zero(bs)
% bs = nan2zero(bs)
%
% Identify regions in obj.parcel_means with NaNs
% Replace all values in identified regions with zero, region-wise
% Needed for some operations and graph metrics
% This could be problematic if NaNs only occur at a few time points
%
% This is a 'stub' function right now, but could be extended to include
% voxel-wise data if available, and have different options for how to
% remove/replace (e.g., element-wise with col means; row- and col-based imputation

len = length(bs.region_dat);

for i = 1:len
    
    rdat = bs.region_dat{i};
    
    whcols = any(isnan(rdat), 1);
    
    if any(whcols)
        
        rdat(:, whcols) = 0;
        
        bs.region_dat{i} = rdat;
        
    end
    
end

end
