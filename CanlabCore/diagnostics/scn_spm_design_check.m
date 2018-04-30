function out = scn_spm_design_check(spm_results_dir, varargin)
% Run in a single-subject (first-level) SPM directory to check 
% design matrix variance inflation factors and high-pass filtering.
% Prints out table of regressors and their above-threshold VIFs (see options).
% Saves .png images of the key figures.
%
% :Usage:
% ::
%
%     scn_spm_design_check(spm_results_dir, varargin)
%
% :Optional Inputs:
%
%   **'events_only':**
%        Show plots and diagnostics for ONLY events, not nuisance covariates or
%        other user-specified regressors.  Useful when you have many nuisance
%        covs.
%
%   **'vif_thresh', t':**
%        Only regressors with a VIF > t will be printed in VIF table.
%
%   **'sort_by_vif'':**
%        Sort regressors in VIF table by VIF (DEFAULT: order regressors as in model).   
%
% Calls: scn_spm_choose_hpfilter.m, scn_spm_get_events_of_interest.m
%
% :Examples:
% ::
%
%    scn_spm_design_check(pwd, 'events_only');
%
% ..
%    Updated: Tor Wager, Aug 2010; Oct 2011: Add 'events_only'; July 2012:
%    fixed for parametric modulators. Luka Ruzic, Sept 2012: added VIF tables.
%    Wani Woo, Apr, 2018: added an output (out) to return vif values 
% ..


if nargin < 1, spm_results_dir = pwd; end

spmfilename = fullfile(spm_results_dir, 'SPM.mat');
if ~exist(spmfilename, 'file')
    error('SPM.mat does not exist in %s\n', spm_results_dir); 
end
load(spmfilename);

VIFTHRESH = 1.3;
EVENTS_ONLY = false;
SORTBYVIF = false;

i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'events_only'
                EVENTS_ONLY = true;
            case 'vif_thresh'
                i=i+1;
                VIFTHRESH = varargin{i};
            case 'sort_by_vif'
                SORTBYVIF = true;
            case 'sort_by_reg'
                % backwards compatibility; do nothing
            otherwise
                error(['Unrecognized argument: ' varargin{i}])
        end
    else
        error(['Unrecognized argument: ' varargin{i}])
    end
    i=i+1;
end


%% Optional inputs


% Gets events of interest: All regressors, or events only if 'events_only'
% is input as keyword
if EVENTS_ONLY
    wh_cols = scn_spm_get_events_of_interest(SPM, 'events_only');
else
    wh_cols = scn_spm_get_events_of_interest(SPM);
end


%%

create_figure('X', 2, 2); 
imagesc(zscore(SPM.xX.X)); set(gca, 'YDir', 'Reverse');
colormap gray
title('Full Design Matrix (zscored)');
axis tight
drawnow

subplot(2, 2, 2)
imagesc(zscore(SPM.xX.X(:, wh_cols))); set(gca, 'YDir', 'Reverse');
colormap gray
title('Design Matrix: Of-interest (zscored)');
axis tight
drawnow

% Variance Inflation Factors for regs of interest (iC)

warning off % some nuisance covs could be redundant; we don't care.
allvifs = getvif(SPM.xX.X(:, SPM.xX.iC), 0);
warning on

subplot(2, 2, 3);
allvifs = allvifs(wh_cols);

plot(allvifs, 'ko', 'MarkerFaceColor', [1 .5 0]);

ylabel('Variance inflation factor'); xlabel('Predictor number');
plot_horizontal_line(1, 'k');
plot_horizontal_line(2, 'b--');
plot_horizontal_line(4, 'r--');
plot_horizontal_line(8, 'r-');
disp('Variance inflation: 1 (black line) = minimum possible (best)');
disp('Successive lines indicate doublings of variance inflation factor.');
title('Var. Inflation (VIFs) in full design');

% Variance Inflation Factors for ONLY of-interest
% now we care if nuisance covs are redundant.
vifs = getvif(SPM.xX.X(:, wh_cols), 0);
subplot(2, 2, 4);
plot(vifs, 'ko', 'MarkerFaceColor', [1 .5 0]);

ylabel('Variance inflation factor'); xlabel('Predictor number');
plot_horizontal_line(1, 'k');
plot_horizontal_line(2, 'b--');
plot_horizontal_line(4, 'r--');
plot_horizontal_line(8, 'r-');
title('VIFs for ONLY of-interest regs');

% table
if SORTBYVIF
    [ignore ord] = sort(allvifs,'descend'); %#ok
else
    ord = 1:numel(allvifs);
end
fprintf('\nVIFs greater than %d\n',VIFTHRESH)
fprintf('\n%7s   %s\n','VIF','REGRESSOR')
for i=1:numel(ord)
    if allvifs(ord(i)) >= VIFTHRESH
        fprintf('%7.2f   reg%04d: "%s"\n',allvifs(ord(i)),wh_cols(ord(i)),SPM.xX.name{wh_cols(ord(i))})
    end
end

try 
    scn_export_papersetup; 
catch
end
saveas(gcf, 'Variance_Inflation', 'png');
disp('Saved Variance_Inflation.png in SPM directory');

if EVENTS_ONLY
    scn_spm_choose_hpfilter(spm_results_dir, 'events_only');
else
    scn_spm_choose_hpfilter(spm_results_dir);
end

%spm_efficiency('SPM.mat');
% saveas(gcf, 'SPM_efficiency', 'png');
% disp('Saved SPM_efficiency.png in current directory');

try 
    scn_export_papersetup(500); 
catch
end
saveas(gcf, 'High_pass_filter_analysis', 'png');
disp('Saved High_pass_filter_analysis.png in SPM directory');

out.allvifs = allvifs;
out.name = SPM.xX.name(wh_cols);

end
