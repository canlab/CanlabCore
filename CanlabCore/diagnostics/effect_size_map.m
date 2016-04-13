% This SCRIPT does an analysis of between- and within-person effect sizes
% in an SPM directory and reports on the effect sizes across voxels, the
% optimal balance of within-person and between-person sample sizes, and
% power for a replication.  It has yet to be developed into a function or
% toolbox, but was used to create several example power analyses in Wager
% and Lindquist's teaching materials.
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2010 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
%
% :Examples:
% ::
%
% %% SETUP Stuff
% 
% % Define a mask from one contrast
% % This contrast is [Nback - Baseline], where N-back is balanced across
% % 2-back and 3-back
% 
% % p < .0001 in localizer
% my_p_img = 'rob_p_0001.img';  % This could be any p-value image; or when using SPM, could be t image
% iimg_threshold(my_p_img, 'thr', [0 .0001], 'outnames', 'rob_p_0001_thresh_0001.img');
% mask = which('rob_p_0001_thresh_0001.img');
% 
% %  mask with brain mask as well
% mask_image(mask, which('brainmask.nii'), mask, 'minmask', .2);
% cl = mask2clusters(mask);
% cluster_orthviews(cl, {[0 1 0]}, 'solid', 'overlay', overlay);
% 
% % Calculate effect size/power from another contrast
% % This is the [3-back - 2-back] contrast
% %cd ../robust0002
% groupcon = which('rob_beta_0001.img');
% groupt = which('rob_tmap_0001.img');
% 
% % Now we need individual-subject MSE images (residual variance images) for
% % the test contrast
% 
% cd('/Volumes/Cort/Imaging/CORT/Analysis-NBACK/Errors_excluded_model2')
% resms_images = filenames('Cort*/ResMS.img', 'absolute', 'char');
% 
% wh_contrast = 2;
% 
% % load sample design matrix and pass in SPM
%cd('/Volumes/Cort/Imaging/CORT/Analysis-NBACK/Errors_excluded_model2/group_robust_triervscontrol_cov/robust0001/')
% go to directory you want to save results in...
%


%%

outputdir = 'Effect_Size_Maps'; 
if ~exist(outputdir, 'dir'), mkdir(outputdir); end
logfile = fullfile(pwd, outputdir, 'Effect_Size_Output.txt');

disp(['Saving output in ' logfile])
diary(logfile);
cd(outputdir);

[volInfo, maskdat] = iimg_read_img(mask, 2);

fprintf('Mask created from: %s\n', volInfo.fname);
fprintf('Mask area: %3.0f voxels\n', volInfo.n_inmask);


con = iimg_get_data(mask, groupcon);
t = iimg_get_data(mask, groupt);

N = size(resms_images, 1);

fprintf('Subjects: N = %3.0f \n', N);

resms = iimg_get_data(mask, resms_images);

% calculate sigma^2 within, assuming all individuals' sigmas are the same
sig2r = mean(resms);  % this is the residual sigma^2, only part of the within-subjects variance

X = SPM.xX.xKXs.X;
[nt, num_cols] = size(X);  % number of within-subject time points and columns

c = SPM.xCon(wh_contrast).c;
desvar = c' * inv(X' * X) * c;  % this is good for calculating observed power, but is dependent on the scaling of the contrast vector and the number of images, nt
sig2w = sig2r .* desvar;  % within-subjects contribution to group cope variance, per subject, factoring in nt

% calculate total effect size, d.  This is signed (pos or neg)
d = t ./ sqrt(N);

% Total variance, sum of between and within
sig2 = con.^2 .* N ./ t.^2;   % variance

% Between-subjects variance estimate 
sig2b = sig2 - sig2w;
sig2b(sig2b < 0) = .0001; % var must be positive...constrain

% Z/t-score measure of effect size - for power calc
% Z is not used here, but appears later with varying number of samples and
% sample size
% Z = con .* sqrt(N) ./ ((sig2b + sig2w) .^ .5);

unit_desvar = desvar .* nt;  % free of nt, per-image design contribution

unit_desvar_norm = desvar .* nt ./ (c'*c); % free of scale of c and nt, but need to divide con by sqrt(c'*c) to use this
con_norm = con ./ sqrt((c'*c));  
sig2w_norm = sig2r .* unit_desvar_norm ./ nt;  % within contribution, c'*c smaller
sig2_norm = con_norm.^2 .* N ./ t.^2;   % variance, (c'*c) smaller than sig2
sig2b_norm = sig2_norm - sig2w_norm;  % between contribution, c'*c smaller than original

fprintf('\nPooled average effect size (d): %3.2f \n', mean(d));
fprintf('\nPooled average residual variance (sig2r): %3.3f \n', mean(sig2r));
fprintf('\nDesign-related variance inflation: %3.3f \n', desvar);
fprintf('\nUnit design-related  variance inflation per observation, sig2d, so that sig2w = sig2r * sig2d / nt: %3.3f \n', unit_desvar);
fprintf('\nUnit design-related variance inflation per observation, normed by c''*c: %3.3f \n', unit_desvar_norm);
fprintf('\nWithin-subjects time points from example SPM.mat: %3.3f \n', nt);
fprintf('\nAverage within-subjects variance contribution to COPE, including design ineff. (sig2w): %3.3f \n', mean(sig2w));
fprintf('\nPooled average unit within-subjects variance per observation (unit_sig2w): %3.3f \n', mean(sig2r * unit_desvar));
fprintf('\nPooled average between-subjects variance (sig2b): %3.3f \n', mean(sig2b));

fprintf('\nThe numbers below for sigma_w and sigma_b can be compared to the calculated map values below')
fprintf('\nCall these stdw and stdb, so that se(COPE)=sqrt(sigw^2/nt + sigb^2 / N) ');
fprintf('\nAverage within-subjects standard deviation per observation, stdw : %3.3f \n', mean(sqrt(sig2r * unit_desvar)));
fprintf('Average between-subjects standard deviation: %3.3f \n', mean(sqrt(sig2b)));


%%
% write maps for key things
iimg_reconstruct_vols(sig2r' .^ .5, volInfo, 'outname', 'residual_std.img');
iimg_reconstruct_vols(sig2w' .^ .5, volInfo, 'outname', 'within_subjects_std.img');
iimg_reconstruct_vols(sig2b' .^ .5, volInfo, 'outname', 'between_subjects_std.img');

% show them
underlay = which('avg152T1.nii');
spm_check_registration(char(underlay, underlay, underlay, underlay));
cl = mask2clusters( 'residual_std.img');
cluster_orthviews(cl, 'add', 'handle', 1);
cl = mask2clusters( 'within_subjects_std.img');
cluster_orthviews(cl, 'add', 'handle', 2);
cl = mask2clusters( 'between_subjects_std.img');
cluster_orthviews(cl, 'add', 'handle', 3);

% key equations

% d = t / sqrt(N)             % effect size, Cohen's d
% sig2 = c^2 * N ./ t^2   % variance

% t or Z = con * sqrt(N) ./ sqrt(sig2b + sig2w)

%% Now tradeoff

% Figure out how sig2w would change given different n within-subjects
% (diff. number of time points)

% overall sig2w = diag( sig2r * c(X'*X)c' )
% sig2w shrinks with k (columns in X) / nt, num. time points
% sig2w = uk/n, where u is unit error var per image
% u = n*sig2w/k
% sig2w_nti = uk/nti, where nti is a new proposed number of time points, i
% assuming k does not change, 
% sig2w varies with 1/n, so u = n*sig2w
% conceptual check: variance scales with 1/n, ste scales with 1/sqrt(n)


% checked: this is the same as sig2r * unit_desvar
unit_sig2w = nt * sig2w;  % unit variance per scan, contribution of measurement noise with one effective image...

% scan hours we have to work with
hours = 60;
hperslot = 1.5;
TR = 2;             % could be "effective TR" 

Ni = 3:100;  % number of subjects
hours_per_subject = hours ./ Ni;

sessions_per_subject = hours_per_subject ./ hperslot;
sessions_per_subject = ceil(sessions_per_subject);

% correct hours because we're losing about 30 mins of setup time and
% structural time per session, or 15 mins for repeat sessions

adj_hours_per_subj = hours_per_subject;
adj_hours_per_subj = adj_hours_per_subj - 30/60;  % lost time in first session
adj_hours_per_subj = adj_hours_per_subj - (15/60 * (sessions_per_subject - 1));

images_per_subject = adj_hours_per_subj .* 3600 ./ TR;  % this is nti above

% correct for effective df
% this is a rough-and-ready estimate based on the SET data (Wager et al.,
% 2009, Part 1)
% it assumes a TR of 2 sec, and an AR(2) model with parameters: [.8711
% -.0334], 215 images translate to 180.5 df on average
df_shrink_factor = .8398;
images_per_subject = images_per_subject .* df_shrink_factor;

% eliminate Ni where we can't actually collect functional data at all on
% the group because we don't have enough time
min_images = num_cols;  % one effective independent data point per column in the design matrix is the minimum

wh = images_per_subject < min_images; 
Ni(wh) = [];
adj_hours_per_subj(wh) = [];
hours_per_subject(wh) = [];
sessions_per_subject(wh) = [];
images_per_subject(wh) = [];

% Ni: candidate number of subjects
% Zi: effect size maps for voxels x Ni
[power_uncorrected05, Zi, power_uncorrected001, d, power_corrected] = deal( zeros(length(con), length(Ni)) );

% Uncorrected power critical t-values, two-tailed
% u_unc_05 = tinv(1 - .05/2, Ni);
% u_unc_001 = tinv(1 - .001/2, Ni);

for i = 1:length(Ni)
    
    sig2wi = unit_sig2w ./ images_per_subject(i);
    
    % Zi: t-values, effect size map (absolute value; no sign info; expected t-value)

    [power_uncorrected05(:, i), Zi(:, i), d(:, i)] = power_from_variance(abs(con)', Ni(i), sig2b', sig2wi', .05);
    
    power_uncorrected001(:, i) = power_from_variance(abs(con)', Ni(i), sig2b', sig2wi', .001);
    
end

% This is an effect size measure, really a t-score, averaged across voxels
Zmean = mean(Zi);

% covert to power estimate, given multiple comparisons
% ---------------------------------------------------------------
% first, estimate smoothness
[fwhm,VRpv] = spm_est_smoothness(spm_vol(resms_images), spm_vol(mask));
fprintf('\nSmoothness (FHWM): %3.1f %3.1f %3.1f \n', fwhm);

rpv = iimg_get_data(mask, 'RPV.img');
rpv(isnan(rpv)) = [];
fprintf('\nRESEL count: %3.1f \n', sum(rpv));

% Corrected threshold (u): min of SPM GRF or Bonf
% --------------------------------------------------------------
clear u ub
for i = 1:length(Ni)
    u(i) = spm_uc_RF(.05, [1 Ni(i) - 1], 'T', sum(rpv), 1);
    ub(i) = spm_uc_Bonf(.05, [1 Ni(i) - 1], 'T', volInfo.n_inmask, 1);
    
end
u = min(u, ub);

% Effective number of independent spatial comparisons (NISC)
% Number of independent tests you'd have to Bonferroni correct
% based on in order to get the corrected p-value threshold
% --------------------------------------------------------------
for i = 1:length(Ni)
    pval_thresh = 1 - tcdf(u, Ni(i) - 1);
end
% .05 / NISC = Bonf based on NISC indep. spatial comparisons
NISC = .05 ./ pval_thresh; 
NISC = mean(NISC(Ni > 20));  % Resels perhaps not corrected for variance due to limited df?
% This is arbitrary, but seems reasonable

%

for i = 1:length(Ni)

    power_corrected(:, i) = 1 - tcdf(u(i) - Zi(:, i), Ni(i) - 1);
       
end




%% at best allocation, save power map

wh = find(mean(power_corrected) == max(mean(power_corrected))); wh = wh(1);

disp('_____________________________________________________')
fprintf('Best balance of subjects and time for %3.0f total scan hours\n', hours);
disp('_____________________________________________________')
fprintf('Subjects: %3.0f \n', Ni(wh));
fprintf('Sessions per subject (max session time of %3.1f hours): %3.0f sessions\n', hperslot, sessions_per_subject(wh));
fprintf('Total scan hours per subject: %3.1f; Functional time: %3.1f hours, or %3.0f min\n', hours_per_subject(wh), adj_hours_per_subj(wh), round(adj_hours_per_subj(wh) .* 60));
 
power_map = power_corrected(:, wh); %1 - tcdf(u(wh) - Zi(:, wh), Ni(wh) - 1);

fprintf('Variance summary for in-mask voxels:\n');
fprintf('Percentile:  \t25th\t50th\t75th\t\n');
fprintf('Within-subj: \t%3.2f\t%3.2f\t%3.2f\t\n', prctile(sig2w, [25 50 75]));
fprintf('Between-subj: \t%3.2f\t%3.2f\t%3.2f\t\n', prctile(sig2b, [25 50 75]));

fprintf('Power for FWE-corrected search:\n');
fprintf('Minimum within search area: %3.0f%%\n', 100 .* min(power_map));
fprintf('Maximum within search area: %3.0f%%\n', 100 .* max(power_map));
fprintf('Mean within search area: %3.0f%%\n', 100 .* mean(power_map));

fprintf('Voxels with 80%% power: %3.0f\n', sum(power_map >= .8));

diary off

%
pmapname = sprintf('power_map_%03dhours_N%03d_ftime%03dmin.img', round(hours), Ni(wh), round(adj_hours_per_subj(wh) .* 60));
iimg_reconstruct_vols(power_map, volInfo, 'outname', pmapname);
cl = mask2clusters(pmapname);
cluster_orthviews(cl, 'add', 'handle', 4);

spm_orthviews_name_axis('Residual std', 1);
spm_orthviews_name_axis('Within-Ss std', 2);
spm_orthviews_name_axis('Between-Ss std', 3);
spm_orthviews_name_axis('FWE corrected power', 4);

%% Power allocation plot

% optimal number
create_figure('Optimal allocation', 2, 2);

plot(Ni, sessions_per_subject, 'k',  'LineWidth', 4);
plot(Ni, adj_hours_per_subj, 'k:',  'LineWidth', 4);
legend({'Sessions per subject' 'Functional scan hours per subject'})
set(gca, 'FontSize', 24);
axis tight

sess_change = diff(sessions_per_subject);
sess_change(sessions_per_subject(2:end) > 2) = 0;  % don't mark for more than 3
sess_change = find(sess_change) + 2;

plot_vertical_line(sess_change);
xlabel('Number of subjects');
title(['Allocation of ' num2str(hours) ' scan hours'])

subplot(2, 2, 3);
plot(Ni, mean(power_corrected), 'k', 'LineWidth', 4);
plot(Ni, mean(power_uncorrected001), '-', 'Color', [.5 .5 .5], 'LineWidth', 4);
plot(Ni, mean(power_uncorrected05), 'k', 'Color', [.75 .75 .75], 'LineWidth', 4);
set(gca, 'FontSize', 24);
xlabel('Number of subjects');
ylabel('Power');
legend({'FWE Corrected' ' p< .001' 'p < .05'})
axis tight

plot_vertical_line(sess_change);
hh = plot_horizontal_line(.8); set(hh, 'LineStyle', ':')

% pie graph of between/within
subplot(2, 2, 2);
sig2wi = unit_sig2w ./ images_per_subject(wh);
set(gca, 'Fontsize', 32)
hh = pie([mean(sig2wi) mean(sig2b)],{'Within','Between'});
colormap gray, axis off, axis image
set(findobj(hh, 'Type', 'Text'), 'FontSize', 32)
hh2 = findobj(hh, 'Type', 'Patch');
set(hh2(1), 'FaceColor', [.3 .3 .3]); 
set(hh2(2), 'FaceColor', [.8 .8 .8]); 

subplot(2, 2, 4);
axis off

try
scn_export_papersetup(600);
saveas(gcf, 'Optimal_allocation_of_hours', 'png');
saveas(gcf, 'Optimal_allocation_of_hours', 'fig');
catch
    disp('Cannot save figure')
end

%% BELOW IS A WHOLE NEW SET OF NORMATIVE CALCULATIONS, NOT BASED ON ANY
% PARTICULAR DATASET; these are in chapter as maps

%% power curves plot - uncorrected -- based on observed

% Set up images and N
% scan hours we have to work with
hours = 60;
hperslot = 1.5;
TR = 2;             % could be "effective TR" 

Ni = 3:100;  % number of subjects
hours_per_subject = hours ./ Ni;

sessions_per_subject = hours_per_subject ./ hperslot;
sessions_per_subject = ceil(sessions_per_subject);

% correct hours because we're losing about 30 mins of setup time and
% structural time per session, or 15 mins for repeat sessions

adj_hours_per_subj = hours_per_subject;
adj_hours_per_subj = adj_hours_per_subj - 30/60;  % lost time in first session
adj_hours_per_subj = adj_hours_per_subj - (15/60 * (sessions_per_subject - 1));

images_per_subject = adj_hours_per_subj .* 3600 ./ TR;  % this is nti above

% correct for effective df
% this is a rough-and-ready estimate based on the SET data (Wager et al.,
% 2009, Part 1)
% it assumes a TR of 2 sec, and an AR(2) model with parameters: [.8711
% -.0334], 215 images translate to 180.5 df on average
df_shrink_factor = .8398;
images_per_subject = images_per_subject .* df_shrink_factor;

% eliminate Ni where we can't actually collect functional data at all on
% the group because we don't have enough time
min_images = 30;  % one effective independent data point per column in the design matrix is the minimum

wh = images_per_subject < min_images; 
Ni(wh) = [];
adj_hours_per_subj(wh) = [];
hours_per_subject(wh) = [];
sessions_per_subject(wh) = [];
images_per_subject(wh) = [];

%
clear power_uncorrected001

within_std = [12:1:72];  % sigma-residual * sigma-design, per observation (image volume)
between_std = [1:.1:10];

conval = 1;
con = repmat(conval, 1, length(within_std));

% each element of 3rd dim is a new between
% each col is sample size
% each row is a new within
[power_uncorrected001, d] = deal(zeros(length(within_std), length(Ni), length(between_std)));

for b = 1:length(between_std)

    for i = 1:length(images_per_subject)

        sig2b = repmat(between_std(b), 1, length(within_std)) .^ 2;
        sig2wi = (within_std .^ 2) ./ images_per_subject(i);

        % Zi: t-values, effect size map (absolute value; no sign info; expected t-value)

        %[power_uncorrected05(:, i), Zi(:, i), d(:, i)] = power_from_variance(abs(con)', Ni(i), sig2b', sig2wi', .05);

        [power_uncorrected001(:, i, b),tmp, d(:, i, b)]  = power_from_variance(abs(con)', Ni(i), sig2b', sig2wi', .001);

    end

end

%create_figure('test', 1, 2); imagesc(power_uncorrected001(:, :, 1)); colorbar; subplot(1, 2, 2); imagesc(power_uncorrected001(:, :, end)); colorbar
%
%%
% Line plot
% -----------------------
cv = 0:.07:.7;
create_figure('plot', 1, 3); 

wi_index = [1 (round(length(within_std) ./ 2)) length(within_std)];

for wi = 1:3

subplot(1, 3, wi)
for b = 1:10
plot(Ni, squeeze(power_uncorrected001(wi_index(wi), :, b)), 'Color', [cv(b) cv(b) cv(b)], 'LineWidth', 2);
end

hours_per_subject = images_per_subject * TR ./ 3600;

set(gca, 'FontSize', 18);
axis tight

sess_change = diff(sessions_per_subject);
sess_change(sessions_per_subject(2:end) > 2) = 0;  % don't mark for more than 3
sess_change = find(sess_change) + 2;

plot_vertical_line(sess_change);

if wi == 3
    for i = 1:10, legstr{i} = [texlabel('sigma_B =') sprintf('%3.1f', between_std(i))];
    end
    legend(legstr, 7);
end

if wi == 2
    str = sprintf('Allocation of %3.0f scan hours\n%s = %3.0f', hours, texlabel('sigma_w'), within_std(wi_index(wi)));
    xlabel('Number of subjects');

else
    str = sprintf('%s = %3.0f', texlabel('sigma_w'), within_std(wi_index(wi)));
end

title(str)
if wi == 1, ylabel('Power'); end
end

disp('Subjects'); disp(Ni(10:10:end))
disp('Hours'); disp(floor(hours_per_subject(10:10:end)))
disp('Min'); disp(round(rem(hours_per_subject(10:10:end), 1) .* 60))
%%

% Contour plot
% -----------------------
create_figure('contour',1, 2)
[m, i] = max(power_uncorrected001, [], 2);
i = squeeze(i);
m = squeeze(m);

[X,Y] = meshgrid(between_std, within_std);
%surf_plot_tor(Ni(i), X, Y, texlabel('sigma_w'), texlabel('sigma_b'), 'Optimal N');

opt_N = Ni(i);
[C, H] = contourf(X, Y, Ni(i));
axis tight
set(gca, 'FontSize', 32, 'XLim', [1 6]);

hh = clabel(C, H);
set(hh, 'FontWeight', 'bold', 'FontSize', 24)

cm = colormap(gray);
%cm = sortrows(cm, -1);
colormap(cm(15:end, :))
xlabel(texlabel('Between-subjects error: sigma_B'));
ylabel(texlabel('Within-subjects error: sigma_w'));
title(sprintf('Sample size (N)'))

subplot(1, 2, 2)
opt_min = round(hours_per_subject(i) .* 60);
[C, H] = contourf(X, Y, opt_min);
axis tight
set(gca, 'FontSize', 32, 'XLim', [1 6]);

hh = clabel(C, H);
set(hh, 'FontWeight', 'bold', 'FontSize', 24)
title(sprintf('Min of scanning/subject'))
xlabel(texlabel('Between-subjects error: sigma_B'));

