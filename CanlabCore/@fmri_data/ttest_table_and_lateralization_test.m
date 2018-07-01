function [roi_table, subj_dat] = ttest_table_and_lateralization_test(imgs, varargin)
% Extracts data from 17 (16) Yeo et al. networks divided into left/right
% hemispheres. Tests the average response in each network (averaged over
% voxels), and tests lateralization (R - L hem) for each. Makes a wedge
% plot and prints a table of results and a short report.
% Currently hard-coded for the 17 networks mapped to Schaefer/Yeo et al. 2018
% Cer Ctx, but could be extended to any symmetric parcellation, and to
% subcortical parcellations.
%
% :Usage:
% ::
%
%     [roi_table, subj_dat] = ttest_table_and_lateralization_test(imgs, ['nomontage'])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2018 Tor Wager
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
% :Inputs:
%
%   **imgs:**
%        an fmri_data object
%
% :Optional Inputs:
%   **none yet:**
%
% :Outputs:
%
%   **roi_table:**
%        Table of mean image values and statistics (matlab table object)
%
%   **subj_dat:**
%        Images x Networks matrix of mean image values within-images
%
% :Examples:
% ::
%
%   imgs = load_image_set('emotionreg');
%   [roi_table, subj_dat] = ttest_table_and_lateralization_test(imgs);
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% Defaults values and optional inputs
% --------------------------------------------------------

domontage = true;
% initalize optional variables to default values here.

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'nomontage', domontage = false;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Display helper functions: Called by later scripts
% --------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);


% Create_wedge_plot and extract subject scores for each network
% ------------------------------------------------------------------------

% Custom colors, mirroring clusters across L/R hem networks:

%[colors1, colors2] = deal(scn_standard_colors(16));
[colors1, colors2] = deal(define_yeo_colors);
colors1 = colors1(1:16);  % note: only 16 regions now...fix?
colors2 = colors2(1:16);

colors = {};
indx = 1;
for i = 1:length(colors1)
    colors{indx} = colors1{i}; colors{indx + 1} = colors2{i}; indx = indx + 2;
end

if domontage
    
    [hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(imgs, 'atlases', {'yeo17networks'}, 'montage', 'colorband_colors', colors);
    
else
    
    [hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(imgs, 'atlases', {'yeo17networks'}, 'colorband_colors', colors);
    
end

labels = labels{1}';                               % network labels
subj_dat = double(output_values_by_region{1});    % subject x network scores

% Note:
% cl = extract_roi_averages(imgs, atlas_obj{1});
% accomplishes the same task as apply_parcellation, returns slightly different values due to interpolation.

%% Lateralization test
% ------------------------------------------------------------------------

% because these are ordered L then R, we can subtract them to get a
% lateralization score.  score so we have R - L
% positive scores are R, negative scores are L

k = size(subj_dat, 2);

wh_left = 1:2:k;
wh_right = 2:2:k;

% Data

L_hem_dat = subj_dat(:, wh_left);
R_hem_dat = subj_dat(:, wh_right);

RL_lat_score = R_hem_dat - L_hem_dat;

% Means and t-ttests
num_tests = size(L_hem_dat, 2) + size(R_hem_dat, 2);
num_tests_lat = size(RL_lat_score, 2);

[R, Tr, Pr, sig_r, bonf_r] = run_ttests(R_hem_dat, num_tests);
[L, Tl, Pl, sig_l, bonf_l] = run_ttests(L_hem_dat, num_tests);
[RvsL, T_rl, P_rl, sig_rl, bonf_rl] = run_ttests(RL_lat_score, num_tests_lat);

unc_r = ~cellfun(@isempty, strfind(sig_r, '*'));
unc_l = ~cellfun(@isempty, strfind(sig_l, '*'));
unc_rl = ~cellfun(@isempty, strfind(sig_rl, '*'));

net_labels = labels(1:2:end);
net_labels = strrep(net_labels, 'LH ', '');

roi_table = table(R, Tr, Pr, sig_r, L, Tl, Pl, sig_l, RvsL, T_rl, P_rl, sig_rl);

%roi_table = table(net_labels, R, L, RvsL, RvsLste, t, p, stars_by_condition);

roi_table.Properties.RowNames = net_labels;

roi_table.Properties.Description = '17 networks, R and L hemispheres, from Schaefer 2018';

disp(roi_table)

%% Summary of results

printhdr('Summary of results');

fprintf('Cont = Control, Default = Default Mode, DorsAttn = Dorsal Attention, \nSalVentAttn = Salience/Ventral Attention, SomMot = Somatomotor, TempPar = Temporal/parietal, \nVisCent = Central visual, VisPeri = Peripheral visual.\n');

fprintf('Two-tailed tests, **** = R/L Bonferroni corrected across %3.0f tests (p < %3.4f), \nlateralization corrected across %3.0f tests (p < %3.4f);\n*** = P < .001; ** = P < .01; * = P < .05; + = P < .10.\n', num_tests, .05 / num_tests,  num_tests_lat, .05 / num_tests_lat);

fprintf('\nBonferroni-corrected results:\n______________________________________\n')

wh = bonf_r & bonf_l & R > 0 & L > 0;
laterality_str = 'bilateral';
direction_str = 'increases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = bonf_r & R > 0;
laterality_str = 'right hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'increases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = bonf_l & L > 0;
laterality_str = 'left hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'increases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

disp(' ')

wh = bonf_r & bonf_l & R < 0 & L < 0;
laterality_str = 'bilateral';
direction_str = 'decreases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = bonf_r & R < 0;
laterality_str = 'right hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'decreases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = bonf_l & L < 0;
laterality_str = 'left hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'decreases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

disp(' ')
lat_n_tests = 16;
fprintf('Lateralization tests corrected across %3.0f regions showed ', lat_n_tests);

wh = bonf_rl & RvsL > 0;
fprintf('\nright-lateralized activation in ');
print_lat_summary_line(wh, net_labels);

wh = bonf_rl & RvsL < 0;
fprintf('Left-lateralized activation was found in ');
print_lat_summary_line(wh, net_labels);


%%
fprintf('\nAt uncorrected thresholds (p < .05 uncorrected):\n______________________________________\n')

wh = unc_r & unc_l & R > 0 & L > 0;
laterality_str = 'bilateral';
direction_str = 'increases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = unc_r & R > 0;
laterality_str = 'right hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'increases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = unc_l & L > 0;
laterality_str = 'left hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'increases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

disp(' ')

wh = unc_r & unc_l & R < 0 & L < 0;
laterality_str = 'bilateral';
direction_str = 'decreases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = unc_r & R < 0;
laterality_str = 'right hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'decreases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

wh = unc_l & L < 0;
laterality_str = 'left hemisphere';  % right hemisphere, left hem, bilateral
direction_str = 'decreases';  % increases or decreases
print_summary_line(wh, net_labels, direction_str, laterality_str)

disp(' ')
lat_n_tests = 16;
fprintf('Lateralization tests corrected across %3.0f regions showed ', lat_n_tests);

wh = unc_rl & RvsL > 0;
fprintf('\nright-lateralized activation in ');
print_lat_summary_line(wh, net_labels);

wh = unc_rl & RvsL < 0;
fprintf('Left-lateralized activation was found in ');
print_lat_summary_line(wh, net_labels);

fprintf('\n----------------------------------------------------------\n');

end % function


%%
function [meandat, t, p, stars_by_condition, bonf_sig] = run_ttests(dat_matrix, num_tests)

meandat = nanmean(dat_matrix)';

bonf_thr = 0.05 ./ num_tests;

% t-tests on each region across images (often subjects)
[h, p, ci, stat] = ttest(dat_matrix, 0, 'tail', 'both');
t = stat.tstat';
p = p';

bonf_sig = p < bonf_thr;

%RvsLste = ste(RL_lat_score)';

% Stars for each region
for j = 1:length(p)
    
    if bonf_sig(j), mystr = '****';
    elseif p(j) < .0015, mystr = '***';
    elseif p(j) < .015, mystr = '**';
    elseif p(j) < .055, mystr = '*';
    elseif p(j) < .105, mystr = '+';
    else mystr = ''; xadj = 0;
    end
    
    stars_by_condition{j, 1} = mystr;
    
end % loop through regions

end % function




function print_summary_line(wh, net_labels, direction_str, laterality_str)


labels_wh = net_labels(wh);

if isempty(labels_wh)
    fprintf('No networks showed %s %s.\n', laterality_str, direction_str)
    
else
    fprintf('Networks with activity %s included %s ', direction_str, laterality_str)
    
    if length(labels_wh) == 1
        fprintf(' %s.\n', labels_wh{end});
        
    else
        for i = 1:length(labels_wh)-1, fprintf('%s, ', labels_wh{i}); end
        fprintf('and %s.\n', labels_wh{end});
    end
    
end


end



function print_lat_summary_line(wh, net_labels)


labels_wh = net_labels(wh);

if isempty(labels_wh)
    fprintf('no networks.\n')
    
else
    if length(labels_wh) == 1
        fprintf(' %s.\n', labels_wh{end});
        
    else
        for i = 1:length(labels_wh)-1, fprintf('%s, ', labels_wh{i}); end
        fprintf('and %s.\n', labels_wh{end});
    end
    
end


end


function c = define_yeo_colors

% as in original yeo lab paper
c = [
    120    18   134
    255     0     0
    70   130   180
    42   204   164
    74   155    60
    0   118    14
    196    58   250
    255   152   213
    220   248   164
    122   135    50
    119   140   176
    230   148    34
    135    50    74
    12    48   255
    0     0   130
    255   255     0
    205    62    78];

c = c ./ 255;

c = mat2cell(c, ones(size(c, 1), 1), 3);

end
