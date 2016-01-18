function add2mask(mask, x, r,varargin)
% Adds or subtracts spheres around x coordinates to/from existing mask
%
% :Usage:
% ::
%
%     add2mask(mask, x, r,varargin)
%
% :Inputs:
%
%   **mask:**
%        is a string filename
%
%   **x:**
%        is n x 3 list of coordinates
%
%   **r:**
%        is radius
%
%   a 4th argument causes us to SUBTRACT!
%
% :Examples:
% ::
%
%    mask = 'insula_from_part1.img'
%    x = [-29.6 28.6 5.3; 37.0 27.5 3.2; -34.9 8.5 -14.8;38.1 10.6 -14.8];
%    add2mask(mask,x,8);
%


add3 = sphere_mask(mask,x,r,'add3.img'); % get spheres

% add or subtract
if length(varargin) == 0
    disp('ADDING spheres to mask')
    tor_spm_mean_ui(str2mat(mask,'add3.img'),'add3.img')
    Q = spm_imcalc_ui('add3.img','add3.img','i1>0');
else
    disp('SUBTRACTING spheres from mask')
    Q = spm_imcalc_ui(str2mat(mask,'add3.img'),'add3.img','i1>0 & ~i2');
end

disp('Written new mask: add3.img')

%Q = spm_imcalc_ui('add3.img','add3.img','i1>0');
%P =%str2mat('add1.img','add2.img','add3.img','smoothed_insula_LR_thresh.img')
%tor_spm_mean_ui(P,'smoothed_insula2.img')
%[v,V,xyz,insula_gray_clusters] = mask_intersection([],'smoothed_insula2.img','smoothed_insula2.img',which('ICBM_brainonly_1mm_seg1.img'));


x = [-29.6 28.6 5.3; 37.0 27.5 3.2; -34.9 8.5 -14.8;38.1 10.6 -14.8];
