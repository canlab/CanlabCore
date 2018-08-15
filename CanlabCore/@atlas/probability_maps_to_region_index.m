function dat = probability_maps_to_region_index(dat)
% Use dat.probability_maps to rebuild integer vector of index labels (dat.dat)
%
% dat = probability_maps_to_region_index(dat)
%

% Start: dat has one image per region, with probability values
% convert to integer vector

[maxval, condf] = max(double(full(dat.probability_maps)),[], 2);   % double is safer

allempty = all(dat.probability_maps == 0, 2) | isnan(maxval);  % some out-of-mask values may get NaNs

condf(allempty) = 0;

dat.dat = int32(condf);

end


