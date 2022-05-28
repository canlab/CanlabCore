function data_obj = movement_regressors(data_obj)
% Create a set of movement regressors with quadratic/derivative
% transformations, and identify high-movement outlier images.
%
% Power (2012) framewise displacement: Uses the average absolute deviation across 6 movement parameters.
% Power 2012: "Spurious but systematic correlations in functional connectivity MRInetworks arisefrom subject motion"
% "FDi = |Δdix| + |Δdiy| + |Δdiz| + |Δαi| + |Δβi| + |Δγi|, where Δdix = d(i − 1)x − dix, 
% and similarly for the other rigid body parameters [dix diy diz αi βi γi]. 
% Rotational displacements were converted from degrees to millimeters by calculating 
% displace- ment on the surface of a sphere of radius 50 mm, which is approximately 
% the mean distance from the cerebral cortex to the center of the head."  Also: FD 
% computed with respect to the volume ~2 seconds previous, with respiration pseudo-motion 
% filtered out of head motion traces (Power 2019). 
% Power 2019: 
% "Breathing can cause the head to move, but changes in lung volume can also cause a particular 
% kind of artifact called pseudomotion, which manifests as a shift of the brain when the lung 
% expands (Brosch et al., 2002; Durand et al., 2001; Raj et al., 2001)."
%
% Yoni's protocol for OLP:
% spike criteria: FD > .25mm OR abs(zscore(DVARS)) > 3, and the following 4 volumes (following ~2 sec)
% bandpass filter, passband = [.01 .1] Hz

mvmt = importdata(mvmt_file);

% fill out 24 movement-related parameters per run
% ----------------------------------------------------------------
mvmt2 = zscore(mvmt) .^ 2;
mvmt3 = [zeros(1, 6); diff(zscore(mvmt))];
mvmt4 = zscore(mvmt3) .^ 2;

% simple geometric mean across translation and rotation; not absolute movement, which varies by voxel
% ----------------------------------------------------------------
geo_displacement = [0; sum(diff(mvmt) .^ 2, 2) .^ .5];

% adjust outliers: add geo displacement

est_outliers_corr(geo_displacement > 0.25) = true;
est_outliers_uncorr(geo_displacement > 0.25) = true;
