cbnames = {'Left I-IV' 'Right I-IV' 'Left V' 'Right V' 'Left VI' 'Vermis VI' 'Right VI' ...
    'Left Crus I' 'Vermis Crus I' 'Right Crus I' 'Left Crus II' 'Vermis Crus II' 'Right Crus II' ...
    'Left VIIb' 'Vermis VIIb' 'Right VIIb' 'Left VIIIa' 'Vermis VIIIa' 'Right VIIIa' ...
    'Left VIIIb' 'Vermis VIIIb' 'Right VIIIb' 'Left IX' 'Vermis IX' 'Right IX' 'Left X' 'Vermis X' 'Right X'};

cbnums = 201:201+length(cbnames)-1;



% load('/Users/tor/Documents/matlab_code_extras/3DheadUtility/Atlases/LPBA40/LBPA40_spm5_label_clusters.mat')
% for i = 1:length(cl)
% loninames{i} = cl{i}(1).shorttitle;
% loninums(i) = 20 + i;
% end

loniname = '/Users/tor/Documents/matlab_code_extras/3DheadUtility/Atlases/LPBA40/maxprob/lpba40_labelID.txt';
[loninums, loninames] = textread(loniname, '%f%s', 'delimiter', '\t');

mask_image = which('atlas_labels_combined.img');

mask = fmri_data(mask_image, which('brainmask.nii'));
mask.additional_info{1} = cbnames;
mask.additional_info{2} = cbnums;

%% 

cl = region(mask_image, 'unique_mask_values');

%%
colors = scn_standard_colors(91);
cluster_orthviews();
for i = 1:length(cl)
cluster_orthviews(cl(i), colors(i), 'add', 'solid');
end

%% Johansen-Berg 2005 connectivity atlas

colors = {[1 .5 0] [1 0 1] [0 1 1] [1 0 0] [0 0 1] [0 1 0] [1 1 0]};

cluster_orthviews();
for i = 1:7
cluster_orthviews(cl(i), colors(i), 'add', 'solid');
end

%%

% for i = 1:7, mycoord = cl(i).XYZmm(:, round(size(cl(i).XYZmm, 2)/2)); spm_orthviews('reposition', mycoord); mycoord
% end

thalnames{1} = 'Thalamus - Primary Motor cortex';
thalnames{2} = 'Thalamus - Sensory Cortex';
thalnames{3} = 'Thalamus - Occipital cortex';
thalnames{4} = 'Thalamus - Prefrontal cortex';
thalnames{5} = 'Thalamus - Premotor cortex';
thalnames{6} = 'Thalamus - Posterior parietal cortex';
thalnames{7} = 'Thalamus - Temporal cortex';

%%

cl = region(mask_image, 'unique_mask_values');

%%
colors = scn_standard_colors(91);
cluster_orthviews();
for i = 1:length(cl)
cluster_orthviews(cl(i), colors(i), 'add', 'solid');
end

%%
mask = fmri_data(mask_image, which('brainmask.nii'));
mask.additional_info{1} = {'Names' 'Cluster indx in region vector'};
mask.additional_info{2} = names;
mask.additional_info{3} = clindx;
mask.Y_descrip = 'Value in mask image';

%%
for i = 1:length(cl)
    
    clindx(i, 1) = i;
    
    mask.Y(i, 1) = mode(cl(i).val);
    
    if mask.Y(i) < 8
        names{i} = thalnames{i};
        
    elseif mask.Y(i) > 200 & mask.Y(i) < 250
        wh = cbnums == mask.Y(i, 1);
        names{i} = cbnames{wh};
        
    elseif any(loninums == mask.Y(i, 1))  %mask.Y(i) > 7 & mask.Y(i) < length(loninames) + 20
        wh = loninums == mask.Y(i, 1);
        names{i} = loninames{wh};

        
    else
        spm_orthviews('reposition', cl(i).mm_center); cl(i).mm_center
        names{i} = input('Name: ', 's');
    end
    
end

%%

cluster_orthviews_montage(6, 'axial');
cluster_orthviews_montage(6, 'sagittal');
han = makelegend(names(1:7), colors(1:7));
han = makelegend(names(64:77), colors(64:length(names)));
han = makelegend(names(78:end), colors(78:length(names)));

%%

pos = [519    53   567   881];

for i = 8:12:63
    han = makelegend(names(i:i+11), colors(i:i+11));
    set(gcf, 'Position', pos, 'Tag', 'tmp');
end

%% Thal

wh = 1:7;

cluster_orthviews();
for i = wh
cluster_orthviews(cl(i), colors(i), 'add', 'solid');
end
han = makelegend(names(wh), colors(wh));
set(gcf, 'Position', pos, 'Tag', 'tmp');
z = cat(2, cl(wh).XYZmm);
z = [min(z(3, :)) max(z(3,:))];

cluster_orthviews_montage(6, 'axial', [], 'range', z);

%% CB

wh = 64:91;

cluster_orthviews();
for i = wh
cluster_orthviews(cl(i), colors(i), 'add', 'solid');
end
han = makelegend(names(wh), colors(wh));
set(gcf, 'Position', pos, 'Tag', 'tmp');
z = cat(2, cl(wh).XYZmm);
z = [min(z(3, :)) max(z(3,:))];

cluster_orthviews_montage(6, 'axial', [], 'range', z);
56
%% Frontal

wh = 8:23;

cluster_orthviews();
for i = wh
cluster_orthviews(cl(i), colors(i), 'add', 'solid');
end
han = makelegend(names(wh), colors(wh));
set(gcf, 'Position', pos, 'Tag', 'tmp');
z = cat(2, cl(wh).XYZmm);
z = [min(z(3, :)) max(z(3,:))];

cluster_orthviews_montage(6, 'axial', [], 'range', z);

%%
save atlas_labels_combined_info mask cl names


%% Create new, simple mask with main regions

mask2 = mask;

newnames = {'Thalamus' 'Cerebellum' 'Frontal_cortex_cingulate' 'Parietal_cortex' ...
    'Occipital_cortex' 'Temporal_cortex' 'Medial_temporal_lobe' 'Insula' ...
    'Basal_ganglia' 'Brainstem'   };
wh = {[1:7] [64:91 62] [8:23 54:55] [24:31] ...
    [32:39] [40:45 48:49] [46:47 50:51 60:61] [52:53] ...
    [56:59] [63]};
    
allvox = zeros(size(mask2.dat));

for i = 1:length(newnames)
    fprintf('%s: ', newnames{i});
    
    for j = wh{i}
        
        fprintf('%3.0f ', j); 
        
        vox = abs(mask.dat - j) < 1000*eps;
        %mask2.dat(vox) = i;
        
        allvox(vox) = i;
    end
    
    fprintf('\n');
end

mask2.dat = allvox;
cl2 = region(mask2, 'unique_mask_values');

%% Not sure what's wrong with the above - try another way

mask2 = mask;
cl3 = region(mask, 'unique_mask_values');

mask2.dat = zeros(size(mask2.dat));

for i = 1:length(newnames)
    fprintf('%s: ', newnames{i});
    
    m = iimg_clusters2indx(cl3(wh{i}), mask.volInfo);
    
    m = m(mask.volInfo.wh_inmask);
    
    mask2.dat(logical(m)) = i;
end

orthviews(mask2);

mask = mask2;
cl = region(mask, 'unique_mask_values');

mask.Y_names = newnames;
for i = 1:length(cl)
cl(i).shorttitle = newnames{i};
end

[dd ff ee] = fileparts(which('atlas_labels_combined.img'));

fname = fullfile(dd, 'atlas_combined_gross_regions.mat');
save(fname, 'mask', 'cl');






