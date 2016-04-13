function cl = extract_indiv_peak_data(cl,imgs)
% Purpose: to find individually significant regions within each subject
% and save average timecourses for each individual ROI for each subject
%
% :Usage:
% ::
%
%    cl = extract_indiv_peak_data(cl,imgs)
%
% :Inputs:
%
%   **cl:**
%        is a clusters structure with one element per region.  Each element
%        (cluster) is a structure containing coordinates and data.
%
%   **imgs:**
%        is a string matrix with one row per subject, with names of images
%        used to define thresholds.  These may be contrast or t-images from
%        individual subjects
% 
% The method of extraction is defined in cluster_tmask, which is currently
% to use 50% of voxels with the highest values in imgs(subject) for each
% subject, or 100% if 50% returns less than 5 voxels.  
% It's easy in cluster_tmask to use absolute t-thresholds instead.
%
% :Required fields of cl:
%
%   **XYZmm:**
%        3 x k list of k coordinates in mm space
%
%   **raw_data:**
%        time x voxels x subjects matrix of data for this region
%
% Outputs appended to cl structure:
% indiv_timeseries, time x subjects averaged individual ROI timecourses
%
% :.INDIV:, a structure with the following information:
%   - tname, the name of the t-mask (or contrast mask) entered
%   - XYZ, voxel coords for significant voxels for each subject
%   - XYZmm, mm coords for sig voxels for each subject
%   - sigt, logical matrix subjects x voxels for sig voxels
%   - maxt, max map value for each subject
%   - center, average coordinate for sig voxels for each subject
%   - mm_center, the same in mm
% 
% The spatial information may be used to correlate peak voxel location with
% behavior, for example.
%
% ..
%    Tor Wager, 7/2005
% ..


disp('extract_indiv_peak_data: getting individual significance regions from raw_data field')
fprintf(1,'Subject ');

for i = 1:size(imgs,1)
    
    fprintf(1,'. %3.0f ',i);
    cl = cluster_tmask(cl,imgs(i,:),i);
    
end

fprintf(1,'\n')

return
