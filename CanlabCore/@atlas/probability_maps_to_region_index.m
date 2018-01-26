function dat = probability_maps_to_region_index(dat)
% Convert dat.probability_maps to dat.dat integer vector or region index values
%
% dat = probability_maps_to_region_index(dat)
%

% Start: dat has one image per region, with probability values
% convert to integer vector

[~, condf] = max(double(dat.probability_maps),[], 2);   % double is safer

allempty = all(dat.probability_maps == 0, 2);

condf(allempty) = 0;

dat.dat = int32(condf);

end


