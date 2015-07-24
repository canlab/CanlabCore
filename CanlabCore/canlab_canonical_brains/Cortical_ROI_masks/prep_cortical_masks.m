dat = fmri_data(which('atlas_labels_combined.img'));


r = region(dat, 'unique_mask_values');
load('/Users/tor/Documents/matlab_code_canlab/trunk/Atlas_Localization_Tools/atlas_labels_combined_info.mat')

%%
% find which regions index groups we're interested

cluster_orthviews(r([8:9 20:21 54:55]))  % medial

cluster_orthviews(r([16:19]))  % lat OFC

cluster_orthviews(r([10:15]))  % lateral

cluster_orthviews(r([52:53]), 'add')  % insula

ctx = fmri_data(which('grey.nii'));
ctx.dat(ctx.dat < .25) = 0;

%% Insula
% =============================================================

wh_regions = [52:53];

u = unique(dat.dat(dat.dat ~= 0))

z = zeros(size(dat.dat));

for i = 1:length(wh_regions)
z(double(dat.dat) == u(wh_regions(i))) = 1;
end

dat_r = dat;
dat_r.dat = z;

dat_r = preprocess(dat_r, 'smooth', 4);
dat_r = threshold(dat_r, [.1 1.1], 'raw-between');
orthviews(dat_r)

rr = region(dat_r);
rr = region2struct(rr);
cluster_orthviews(rr, 'unique')

dat_r.fullpath = '/Users/tor/Documents/matlab_code_canlab/trunk/Atlas_Localization_Tools/insula_liberal_no_operc.img';
write(dat_r);

% save
dat_insula = dat_r;
dat_insula = replace_empty(dat_insula);
is_insula = dat_insula.dat > 0;

% masked version - cortical ribbon

dat = fmri_data(which('insula_liberal_no_operc.img'));

mask = which('SPM8_colin27T1_seg.img');
mask = fmri_data(mask);
orthviews(mask)
mask = threshold(mask, [105 140], 'raw-between');

dat = apply_mask(dat, mask);
orthviews(dat)


%% Medial
% =============================================================

wh_regions = [8:9 20:21 54:55];

u = unique(dat.dat(dat.dat ~= 0))

z = zeros(size(dat.dat));

for i = 1:length(wh_regions)
z(double(dat.dat) == u(wh_regions(i))) = 1;
end

dat_r = dat;
dat_r.dat = z;


mmcut = mm2voxel([0 -24 0], dat_r.volInfo);
mmcut = mmcut(2);

omit = dat_r.volInfo.xyzlist(:, 2) < mmcut;
dat_r.dat(omit) = 0;

dat_r.dat(~ctx.dat) = 0;

rr = region(dat_r);
rr = region2struct(rr);
cluster_orthviews(rr, 'unique')

dat_r.fullpath = '/Users/tor/Documents/matlab_code_canlab/trunk/Atlas_Localization_Tools/medial_frontal.img';
write(dat_r);

%% Orbital
% =============================================================

wh_regions = [16:19];

u = unique(dat.dat(dat.dat ~= 0))

z = zeros(size(dat.dat));

for i = 1:length(wh_regions)
z(double(dat.dat) == u(wh_regions(i))) = 1;
end

dat_r = dat;
dat_r.dat = z;

rr = region(dat_r);

rr = region2struct(rr);

cluster_orthviews(rr, 'unique')

dat_r.fullpath = '/Users/tor/Documents/matlab_code_canlab/trunk/Atlas_Localization_Tools/lateral_orbitofrontal.img';
write(dat_r);

%% Lateral PFC
% =============================================================

wh_regions = [8:9 10:15];

u = unique(dat.dat(dat.dat ~= 0));

z = zeros(size(dat.dat));

for i = 1:length(wh_regions)
z(double(dat.dat) == u(wh_regions(i))) = 1;
end

dat_r = dat;
dat_r.dat = z;



mmcut = mm2voxel([-12 0 0], dat_r.volInfo);
mmcut = mmcut(1);
b1 = mmcut;

mmcut = mm2voxel([12 0 0], dat_r.volInfo);
mmcut = mmcut(1);
b2 = mmcut;

omit = dat_r.volInfo.xyzlist(:, 1) < b1 & dat_r.volInfo.xyzlist(:, 1) > b2;
dat_r.dat(omit) = 0;


mmcut = mm2voxel([0 0 -11], dat_r.volInfo);
mmcut = mmcut(3);

omit = dat_r.volInfo.xyzlist(:, 3) < mmcut;
dat_r.dat(omit) = 0;


% omit z < 45 and x < 19
mmcut = mm2voxel([-19 0 45], dat_r.volInfo);
b1 = mmcut;

mmcut = mm2voxel([19 0 45], dat_r.volInfo);
b2 = mmcut;

omit = dat_r.volInfo.xyzlist(:, 3) < b1(3) & dat_r.volInfo.xyzlist(:, 1) < b1(1) & dat_r.volInfo.xyzlist(:, 1) > b2(1);
dat_r.dat(omit) = 0;


% omit y < 14 and x < 19
mmcut = mm2voxel([-19 14 0], dat_r.volInfo);
b1 = mmcut;

mmcut = mm2voxel([19 14 0], dat_r.volInfo);
b2 = mmcut;

omit = dat_r.volInfo.xyzlist(:, 2) < b1(2) & dat_r.volInfo.xyzlist(:, 1) < b1(1) & dat_r.volInfo.xyzlist(:, 1) > b2(1);
dat_r.dat(omit) = 0;

% omit insula
dat_r.dat(is_insula) = 0;

% omit non-cortex

dat_r.dat(~ctx.dat) = 0;

% omit stray voxels
rr = region(dat_r);

rr(cat(1, rr.numVox) < 100) = [];
dat_r = region2imagevec(rr);



rr = region(dat_r);

rr = region2struct(rr);

cluster_orthviews(rr, 'unique')

dat_r.fullpath = '/Users/tor/Documents/matlab_code_canlab/trunk/Atlas_Localization_Tools/lateral_prefrontal.img';
write(dat_r);