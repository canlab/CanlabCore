function NORM_CHECK = scnlab_norm_check(template, wanat_files, mean_func_files, subjects)
% Compares the similarity of one or two sets of images (wanat_files,
% mean_func_files) to a template image and to one another (via Malanobis
% distance) to determine whether some images are potential outliers.
% This is used to check the quality of spatial warping/normalization for a
% group of subjects, though it could be used for other purposes as well.
%
% :Usage:
% ::
%
%     NORM_CHECK = scnlab_norm_check(template, wanat_files, mean_func_files, subjs)
%
% :Inputs:
%
%   **template:**
%        Char array with name of image of normalization template
%
%   **wanat_files:**
%        Warped (to template) anatomical file names
%
%   **mean_func_files:**
%        Names of mean functional images
%        These images should all be in the same space/in register.  
%
%   **Subjs:**
%        Optional cell array of names for each subject, for display
%        purposes
%
% :Outputs:
%
%   A structure with metrics (NORM_CHECK)
%
%   **NORM_CHECK.global_t1:**
%        global values of first image series (wanat_files)
%
%   **NORM_CHECK.std_t1:**
%        spatial standard deviation of first image series (wanat_files)
%
%   **NORM_CHECK.names_t1:**
%        Names for columns of NORM_CHECK.norm_vs_template
%
%   **NORM_CHECK.subjects:**
%        Cell array of names for each subject
%
%   **NORM_CHECK.norm_vs_template:**
%        Similarity data for subjects (rows) x metrics (cols)
%           {'Dist. from group, actual chi2', 'Mutual info with template', 'Correlation with template'};
%
% NB: Leave mean_func_files empty (e.g., []) to only check structural images
%
% Computes metrics on the goodness of normalization based on multivariate distance,
% mutual information, and correlation with template. Automatically saves a .mat file
% of the results into the current directory.
%
% the template file (i.e., avg152T1.nii) must be in the CURRENT working
% directory and have read/write permissions
%
% USES the subfunction compare_subjects, which may be useful as a
% stand-alone function.
%
% USED in canlab_preproc_norm_check.m

parse_inputs();

if(~exist('subjects', 'var') || isempty(subjects))
    subjects = cell(1, length(wanat_files));
    for i = 1:length(subjects)
        subjects{i} = sprintf('Subj %d', i);
    end
end

[ds, g, mystd, d, d2, c, c2, mi] = compare_subjects(char(wanat_files), mask, 0, 'Normed T1 vs. template', 1, subjects, template);

NORM_CHECK.global_t1 = g';
NORM_CHECK.std_t1 = mystd';
NORM_CHECK.names_t1 = {'Dist. from group, actual chi2', 'Mutual info with template', 'Correlation with template'};
NORM_CHECK.subjects = subjects;
NORM_CHECK.norm_vs_template = [ds(:,1) mi' c'];

disp('Structural image quality values saved in NORM_CHECK.norm_vs_template')
disp('Chi2 from group (higher is worse)  Mut info with template(hi=good)  Corr with template(hi=good)')

try
    save('NORM_CHECK', 'NORM_CHECK');
catch
    disp('CANNOT save NORM_CHECK.mat. Permissions problems?');
end

try
    saveas(gcf, 'Norm_vs_Templ.fig');
    saveas(gcf, 'Norm_vs_Templ.png');
    close(gcf);
catch
    disp('Cannot save images');
end

if(exist('mean_func_files', 'var') && ~isempty(mean_func_files))
    mean_func_files = char(mean_func_files);
    
    anybad = check_spm_matfiles(mean_func_files);
    if anybad, error('You must fix these subjects before proceeding.'); end
    
    [ds, g, mystd, d, d2, c, c2, mi] = compare_subjects(char(mean_func_files), mask, 0, 'Mean norm funct vs. group', 1, subjects);
    
    tmp = sortrows(ds, 2);
    tmp = tmp(:, 1) - tmp(:, 3);
    NORM_CHECK.global = g';
    NORM_CHECK.stdev = mystd';
    NORM_CHECK.names = {'Dist. from group, actual chi2', 'Mutual info with grp. mean', 'Correlation with grp. mean'};
    NORM_CHECK.funct_vs_group = [tmp mi' c'];
    
    disp('Functional image quality values saved in NORM_CHECK.funct_vs_group')
    disp('Chi2 from group (higher is worse)  Mut info with template(hi=good)  Corr with template(hi=good)')
    
    
    try
        saveas(gcf, 'subjplot_Func_vs_Group.fig');
        saveas(gcf, 'subjplot_Func_vs_Group.png');
        close(gcf);
    catch
        disp('Cannot save images');
    end
end
try
    save NORM_CHECK -append NORM_CHECK
catch
    disp('CANNOT save NORM_CHECK.mat. Permissions problems?');
end

plot_results(NORM_CHECK, template, wanat_files, mean_func_files, subjects);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inline functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function parse_inputs()
        switch(spm('Ver'))
            case {'SPM5' 'SPM8' 'SPM12'}
                mask = which('avg152T1.nii');
                if isempty(mask)
                    mask = spm_select(1, 'image', 'Cannot find brain mask: please select', [], pwd());
                end
                if(~exist('template', 'var') || isempty(template)) %#ok<NODEF>
                    template = spm_select(1, 'image', 'No template image passed in: please select', [], pwd());
                end
                if(~exist('wanat_files', 'var') || isempty(wanat_files)) %#ok<NODEF>
                    wanat_files = spm_select(Inf, 'image', 'No subjects'' normalized anatomical images passed in: please select', [], pwd(), 'w.*');
                end
            case 'SPM2'
                mask = which('avg152T1.img');
                if isempty(mask)
                    mask = scan_get_files(1, '*img', 'Cannot find brain mask: please select');
                end
                if(~exist('template', 'var') || isempty(template)) %#ok<NODEF>
                    template = spm_get(1, '*.img', 'No template image passed in: please select');
                end
                if(~exist('wanat_files', 'var') || isempty(wanat_files)) %#ok<NODEF>
                    wanat_files = spm_get(Inf, 'w*img', 'No subjects'' normalized anatomical images passed in: please select');
                end
            otherwise
                error('Unknown SPM version: %s\n', spm('Ver'));
        end
        
        if(~exist('mean_func_files', 'var'))
            mean_func_files = [];
        end
        
        wanat_files = cellstr(wanat_files);
        template = char(template);
    end
end




function plot_results(NORM_CHECK, template, wanat_files, mean_func_files, subjects)
tmp = scale(NORM_CHECK.norm_vs_template); tmp2 = tmp;

% save best two and worst two overall
overall = tmp; overall(:, 2:3) = -1*overall(:, 2:3); overall = mean(overall, 2);   % higher is worse
NORM_CHECK.overall_t1 = overall; overall2 = overall;
tmp = find(overall == min(overall)); wh(1) = tmp; overall(tmp) = NaN;
tmp = find(overall == min(overall)); wh(2) = tmp; overall(tmp) = NaN;
tmp = find(overall == max(overall)); wh(3) = tmp; overall(tmp) = NaN;
tmp = find(overall == max(overall)); wh(4) = tmp; overall(tmp) = NaN;
P = strvcat(template, wanat_files{wh});
NORM_CHECK.best_t1 = P(1:2, :);
NORM_CHECK.worst_t1 = P(3:4, :);

try
    spm_check_registration(P);
    
    h = sort(get(gcf, 'Children'));
    axes(h(3)); title('Canonical template');
    axes(h(7)); title(['Best: ' subjects{wh(1)}]);
    axes(h(11)); title(['Best: ' subjects{wh(2)}]);
    axes(h(15)); title(['Worst: ' subjects{wh(3)}]);
    axes(h(19)); title(['Worst: ' subjects{wh(4)}]);
    
    axes('Position', [.6 .05 .3 .25]);
    plot(tmp2, 1:length(tmp2));
    try legend(NORM_CHECK.names_t1, 'Location', 'Best'); catch, end
    hold on; plot(overall2, 1:length(overall), 'ko-', 'LineWidth', 2);drawnow
    set(gca, 'YTickLabel', subjects, 'YTick', 1:length(overall))
    xlabel('Value (higher is worse for black and blue lines)')
catch
end

disp('Line Plot: For Structural vs. template')
disp('Subject order goes from bottom to top')

drawnow
saveas(gcf, 'Best_and_worst_t1.fig');
saveas(gcf, 'Best_and_worst_t1.png');

% Separate line plot
create_figure('Norm_quality_anatomical');
plot(tmp2, 1:length(tmp2));
try legend(NORM_CHECK.names_t1, 'Location', 'Best'); catch, end
hold on; plot(overall2, 1:length(overall), 'ko-', 'LineWidth', 2);drawnow
set(gca, 'YTickLabel', subjects, 'YTick', 1:length(overall))
xlabel('Value (higher is worse for black and blue lines)')
saveas(gcf, 'Norm_quality_anatomical.png');

if(exist('mean_func_files', 'var') && ~isempty(mean_func_files))
    tmp = scale(NORM_CHECK.funct_vs_group);
    tmp2 = tmp;
    
    % save best two and worst two overall
    overall = tmp; overall(:, 2:3) = -1*overall(:, 2:3); overall = mean(overall, 2);   % higher is worse
    NORM_CHECK.overall_funct = overall; overall2 = overall;
    tmp = find(overall == min(overall)); wh(1) = tmp; overall(tmp) = NaN;
    tmp = find(overall == min(overall)); wh(2) = tmp; overall(tmp) = NaN;
    tmp = find(overall == max(overall)); wh(3) = tmp; overall(tmp) = NaN;
    tmp = find(overall == max(overall)); wh(4) = tmp; overall(tmp) = NaN;
    if ischar(mean_func_files)      %Schafer edit, mean_func needs to be a cell for the following code
        szmf = size(mean_func_files);
        for idx = 1:szmf(1)
            tmpfunc{idx,1} = mean_func_files(idx,:);
        end
        mean_func_files = tmpfunc;
    end                             %end Schafer edit 10/18/2010
    
    P = strvcat(template, mean_func_files{wh});
    NORM_CHECK.best_funct = P(1:2, :);
    NORM_CHECK.worst_funct = P(3:4, :);
    
    try
        spm_check_registration(P);
        h = sort(get(gcf, 'Children'));
        axes(h(3)); title('Canonical template');
        axes(h(7)); title(['Best: ' subjects{wh(1)}]);
        axes(h(11)); title(['Best: ' subjects{wh(2)}]);
        axes(h(15)); title(['Worst: ' subjects{wh(3)}]);
        axes(h(19)); title(['Worst: ' subjects{wh(4)}]);
        
        axes('Position', [.6 .05 .3 .25]);
        plot(tmp2, 1:length(tmp2));
        try legend(NORM_CHECK.names_t1, 'Location', 'Best'); catch, end
        hold on; plot(overall2, 1:length(overall), 'ko-', 'LineWidth', 2);drawnow
        set(gca, 'YTickLabel', subjects, 'YTick', 1:length(overall))
        xlabel('Value (higher is worse for black and blue lines)')
    catch
    end
    
    saveas(gcf, 'Best_and_worst_funct.fig');
    saveas(gcf, 'Best_and_worst_funct.png');
    
    % Separate line plot
    create_figure('Norm_quality_mean_functional');
    plot(tmp2, 1:length(tmp2));
    try legend(NORM_CHECK.names_t1, 'Location', 'Best'); catch, end
    hold on; plot(overall2, 1:length(overall), 'ko-', 'LineWidth', 2);drawnow
    set(gca, 'YTickLabel', subjects, 'YTick', 1:length(overall))
    xlabel('Value (higher is worse for black and blue lines)')
    saveas(gcf, 'Norm_quality_mean_functional.png');
    
end % if mean func entered

end % plot results
