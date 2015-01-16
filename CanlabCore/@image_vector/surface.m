function [all_surf_handles, pcl, ncl] = surface(obj)
% [all_surf_handles, pcl, ncl] = surface(obj)
%
% Examples:
% ------------------------------------------------------------------------
% % create an initial surface plot from an fmri_data object:
% han = surface(regionmasks{2});  
%
% Now add a second region in green:
% cluster_surf(region(regionmasks{2}), {[0 1 0]}, han, 5);

if size(obj.dat, 2) > 1
    obj = mean(obj);
end

r = region(obj);


[all_surf_handles, pcl, ncl] = surface(r);

end

