function canlab_create_wm_ventricle_masks(wm_mask, gm_mask, varargin)
% This function saves white matter and ventricle masks.
%
% :Usage:
% ::
%
%     function canlab_create_wm_ventricle_masks(wm_mask, gm_mask)
%
% :Inputs:
%
%   **wm_mask:**
%        white matter structural image file
%        ::
%
%            wm_mask = filenames('Structural/SPGR/wc2*.nii', 'char', 'absolute');
%
%   **gm_mask:**
%        gray matter structural image file 
%        ::
%
%            gm_mask = filenames('Structural/SPGR/wc1*.nii', 'char', 'absolute');
%
% :Optional:  
%    You can specify how liberal or conservative to be in estimating white
%    matter and ventricles. 1 is most conservative and will yield no
%    ventricles, and 0 is very liberal. The default of wm_thr is .9, the
%    default of vent_thr is .9. 
%    e.g)
%      - 'wm_thr', .99
%      - 'vent_thr', .95
%
% :Output:
%
%   **"white_matter.img" and "ventricles.img":**
%        in the same folder of the input structural files
%
% ..
%    5/4/2012 by Tor Wager and Wani Woo
%    7/16/2014 creation of ventricle mask updated by Yoni Ashar
%    12/11/2014 fixed some bugs by Wani Woo
% ..

canonvent_mask = which('canonical_ventricles.img');
bstem = which('spm2_brainstem.img');
canonical_wm = which('white.nii');

if isempty(canonvent_mask) || isempty(bstem) || isempty(canonical_wm)
    error('If you want to use this function, you need ''canonical_ventricles.img'', ''spm2_brainstem.img'', and ''white.nii'' in your path.');
end

if ~exist(wm_mask, 'file') || ~exist(gm_mask, 'file')
    error('Mask files passed in as parameters do not exist.')
end

vent_thr = .1; % 1 - 0.8
wm_thr = .9;
if any(strcmp(varargin, 'wm_thr'))
    wm_thr = varargin{find(strcmp(varargin, 'wm_thr'))+1};
end

if any(strcmp(varargin, 'vent_thr'))
    vent_thr = 1 - varargin{find(strcmp(varargin, 'vent_thr'))+1};
end

%% WHITE MATTER

wm = statistic_image('image_names', wm_mask);
wm = threshold(wm, [.99 1.1], 'raw-between');
%orthviews(wm)

% mask with canonical
canonwm = statistic_image('image_names', canonical_wm);
canonwm = threshold(canonwm, [wm_thr 1.1], 'raw-between');
canonwm.dat(~canonwm.sig) = 0;

wm = apply_mask(wm, canonwm);
%%

bstem = statistic_image('image_names', bstem);
bstem = resample_space(bstem, wm, 'nearest');

wm = replace_empty(wm);
bstem = replace_empty(bstem); 
% remove brainstem voxels
wm = remove_empty(wm, logical(bstem.dat));

% write
d = fileparts(wm.fullpath);
wm.fullpath = fullfile(d, 'white_matter.img');
wm.dat(~wm.sig) = 0;
write(wm);


%% VENTRICLE


% find voxels in canonical brain and in canonical ventricles,
% but not in WM or gray matter

gm = statistic_image('image_names', gm_mask);
gm = threshold(gm, [vent_thr 1.1], 'raw-between'); 

wm = statistic_image('image_names', wm_mask);
wm = threshold(wm, [vent_thr 1.1], 'raw-between'); 

% before reconstructing, "manually" enforce the threshold
gm.dat = gm.dat .* gm.sig;
wm.dat = wm.dat .* wm.sig;

wm_r = reconstruct_image(wm);
gm_r = reconstruct_image(gm);

in_neither = ~(gm_r | wm_r);
vent = rebuild_volinfo_from_dat(image_vector('image_names', wm_mask), in_neither(:)); % wm_mask just provides image space

% now only take part of image in the canonical ventricles area
vent = apply_mask(vent, statistic_image('image_names', canonvent_mask));

% write
vent.fullpath = fullfile(fileparts(wm.fullpath), 'ventricles.img');
write(vent)

return
