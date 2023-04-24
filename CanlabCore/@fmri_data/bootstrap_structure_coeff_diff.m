function [r1_obj, r2_obj, rdiff_obj] = bootstrap_structure_coeff_diff(obj, pattern1, pattern2)
% nps = load_image_set('nps');
% vps = load_image_set('vps');
% obj = load_image_set('emotionreg');
%
% pattern1 = nps;  pattern2 = vps;
% [r1_obj, r2_obj, rdiff_obj] = bootstrap_structure_coeff_diff(obj, pattern1, pattern2)
%
% % Re-threshold maps and display again:
% r1_obj = threshold(r1_obj, .05, 'unc');
% rdiff_obj = threshold(rdiff_obj, .05, 'unc');


obj.dat = double(obj.dat);  % superstition

%% Apply patterns

pexp1 = double(apply_mask(obj, pattern1, 'pattern_expression', 'cosine_similarity'));

pexp2 = double(apply_mask(obj, pattern2, 'pattern_expression', 'cosine_similarity'));

%% Calculate structure coefficient maps for full sample

r1 = corr(obj.dat', pexp1);

r2 = corr(obj.dat', pexp2);

rdiff = r1 - r2;


%% Set up bootstrap samples

nboot = 100;
bootsam = setup_boot_samples(pexp1, nboot);

%%

nvox = size(obj.dat, 1);
[r1_boot, r2_boot, rdiff_boot] = deal(NaN .* zeros(nvox, nboot));

tic

for i = 1:nboot

    indx = bootsam(:, i);

    obj_boot = get_wh_image(obj, indx);

    pexp1_boot = pexp1(indx);
    pexp2_boot = pexp2(indx);

    % Get summary statistic maps for this bootstrap sample
    % voxels x bootstrap samples, each column is a map
    r1_boot(:, i) = corr(double(obj_boot.dat'), pexp1_boot);

    r2_boot(:, i) = corr(double(obj_boot.dat'), pexp2_boot);

    rdiff_boot(:, i) = r1_boot(:, i) - r2_boot(:, i);

end % bootstrap iteration

toc

%%

r1_obj = get_statistic_image(r1_boot, r1, obj, 'Structure coeff map for pattern 1');

r2_obj = get_statistic_image(r2_boot, r2, obj, 'Structure coeff map for pattern 2');

rdiff_obj = get_statistic_image(rdiff_boot, rdiff, obj, 'Structure coeff map for pattern 1 - 2');


end % main function



%% Subfunctions

function r1_obj = get_statistic_image(r1_boot, r1, obj, descrip)

r1_p = bootbca_pval(0, @mean, r1_boot' ,r1', obj.dat');
r1_p = r1_p';
% r1_z = r1_z';

r1_fdr_thresh = FDR(r1_p, 0.05);
if isempty(r1_fdr_thresh),  r1_fdr_thresh = -Inf; end

r1_obj = statistic_image('dat', r1, ...
    'volInfo', obj.volInfo, ...
    'p', r1_p, ...
    'sig', r1_p < r1_fdr_thresh, ...
    'ste', [], ...
    'dat_descrip', descrip, ...
    'removed_voxels', obj.removed_voxels);

end
