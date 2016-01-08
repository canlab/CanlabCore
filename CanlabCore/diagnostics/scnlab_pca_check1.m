function scnlab_pca_check1(imgs, realign_files, X, spersess)
% :Usage:
% ::
%
%     function scnlab_pca_check1(imgs, realign_files or params (t x 6) across all runs, X, spersess)
%
%
% :Inputs:
%
%   **imgs:**
%        list of all image names in order
%
%   **realign_files:**
%        movement param file for each session, names in a cell array, OR
%
%        a t x 6 matrix of realignment parameters across all sessions
%
%   **X:**
%        design matrix; no intercept is needed

% :Examples:
% ::
%
%    % setup code for auditory oddball data
%    cd('/Users/tor/Documents/Tor_Documents/Coursework_and_Teaching/Mind_Res_Net_fMRI_Course_2008/data/auditory_oddball/2subjects-processed/s01/')
%
%    imgs = filenames('*/sw*img','absolute','char')
%    realign_files = filenames('*/rp*txt')
%
%    % LOAD TASK ONSETS and CREATE DESIGN MATRIX
%    onsets{1} = load('novel_stimuli_run1.asc');
%    onsets{2} = load('target_stimuli_run1.asc');
%    onsets{3} = load('standard_stimuli_run1.asc');
%    onsets{4} = load('novel_stimuli_run2.asc');
%    onsets{5} = load('target_stimuli_run2.asc');
%    onsets{6} = load('standard_stimuli_run2.asc');
%
%    regs_per_sess = 3;
%    nsess = 2;
%    for i = 1:length(onsets), onsets{i} = onsets{i}'; end
%    X = cell(1, nsess);
%    X{1} = onsets2delta(onsets(1:3), 1, 249);
%    X{1} = X{1}(:, 1:end-1);
%    X{2} = onsets2delta(onsets(4:6), 1, 249);
%    X{2} = X{2}(:, 1:end-1);
%    X = blkdiag(X{:});

% ..
%    PCA ANALYSIS USING FMRISTAT
%
%    intercept, for later removal
% ..
Xint = intercept_model(spersess, 1:10);

scum = cumsum(spersess);

if ~exist('qc_images', 'dir'), mkdir('qc_images'); end

figure; set(gcf,'Position',[307         310        1109         470]); 
mask_thresh = fmri_mask_thresh(imgs(1,:));
% [V, D] = pca_image(imgs, [], 5, imgs(1,:), mask_thresh, 'pca_raw', 0);

%print('-dpsc2', '-append', 'qc_report');
%scn_export_papersetup(800); saveas(gcf,['qc_images' filesep 'pca_image_raw'],'png');


[V, D] = pca_image(imgs, [], 5, imgs(1,:), mask_thresh, 'pca_noint', 0, Xint);
subplot(1, 2, 1);
for i = 1:length(spersess), plot_vertical_line(scum(i)); end
scn_export_papersetup(800); saveas(gcf,['qc_images' filesep 'pca_image_no_intcpt'],'png');

% -------------------------------------------------------------------------
% LOAD REALIGNMENT PARAMS
% -------------------------------------------------------------------------
if iscell(realign_files)
    disp('LOADING REALIGNMENT PARAMS')

    for i = 1:length(realign_files), rparams{i} = load(realign_files{i}); end;
    rp = cat(1, rparams{:});

    % get RP by session -- for relation to PCs
    for i = 1:length(rparams), rparams{i} = zscore(rparams{i}); end
    rpz = blkdiag(rparams{:});

else
    disp('I THINK MOVEMENT PARAMS ARE ENTERED IN A SINGLE MATRIX (T X 6)')
    rp = realign_files;
    rpz = rp;

end

create_figure('Realignment params'); 
plot(rp);
for i = 1:length(spersess), plot_vertical_line(scum(i)); end

%OPT = scnlab_outlier_id('setup', 'tr', 2, 'spersess', spersess, 'dummy', 1:3, 'hp', 100, 'mad', 5, 'niter', 3, 'mvmt', rpzcat);
OPT = scnlab_outlier_id('setup', 'tr', 2, 'spersess', spersess, 'hp', 100, 'mad', 5, 'niter', 3, 'mvmt', rp);



m = size(V, 2); n = size(rpz, 2); p = size(X, 2);
t = size(V, 1);

% ------------------------------------------------------------------------
% Movement plot
% -------------------------------------------------------------------------
fprintf('%% Var Explained by mvmt:\n');
% for j = 1:4
%     
%     [b, bint, r, rint, stats] = regress(V(:, j), rpz);
%     fprintf('\tComp %3.0f : %3.0f%%\n', j, stats(1)*100);
% end

create_figure('Scatterplots');

[H, Ax,BigAx] = plotmatrix(rp, V);
xlabel('Movement param') ; ylabel('Component')

print('-dpsc2', '-append', 'qc_report');
scn_export_papersetup(800); saveas(gcf,['qc_images' filesep 'motion_component_scatter'],'png');


Xm = intercept(rpz, 'end');
rsquare = zeros(1, 4);

disp('Variance in each component explained by realignment parameters:')

for j = 1:m
    % Center component to avoid counting intercept in r-square?
    % Don't have to, doesn't matter because var operator is 2nd moment
    y = V(:, j);  
    b =  Xm \ y;

    fits = Xm * b;

    rsquare(j) = var(fits) / var(y);
    
    fprintf('\tComp %3.0f : %3.0f%%\n', j, rsquare(j)*100);
end

%%
% -------------------------------------------------------------------------
% Get VIFs
% -------------------------------------------------------------------------

% Task: Variance inflation factors
vifs = getvif(X);

% Variance inflation factors, considering movement params too
vifsm = getvif([X rpz]);
vifsm = vifsm(1:size(X,2));

% -------------------------------------------------------------------------
% PLOT OF RELATIONSHIPS
% -------------------------------------------------------------------------

cc = corrcoef([rpz V]);
mbyc = cc(1:n, n+1:end);

cc = corrcoef([rpz X]);
mbyt = cc(1:n, n+1:end);

cc = corrcoef([X V]);
tbyc = cc(1:p, p+1:end);

create_figure('Relationships', 3, 3);
for i = 1:9, axh(i) = subplot(3, 3, i); end

% % for i = 1:9
% %     pp = get(axh(i), 'Position'); pp(2) = pp(2) - .05; pp(4) = pp(4) + .05; pp(3) = pp(3) + .05; set(axh(i), 'Position', pp);
% % end
axes(axh(1))
axis off

% marginals
axes(axh(2)); cla
plot_matrix_cols(V, 'vertical')
set(gca,'YLim', [0 t+1], 'XLim', [0 m+1], 'XAxisLocation', 'top');
title('Component'), drawnow

axes(axh(3)); cla
plot_matrix_cols(X, 'vertical')
set(gca,'YLim', [0 t+1], 'XLim', [0 p+1], 'XAxisLocation', 'top');
title('Task Design'), drawnow

axes(axh(4)); cla
plot_matrix_cols(rpz)
set(gca,'YLim', [0 n+1], 'XLim', [0 t+1]);
ylabel('Movement'), drawnow

axes(axh(7)); cla
plot_matrix_cols(X)
set(gca,'YLim', [0 p+1], 'XLim', [0 t+1]);
ylabel('Task Design'), drawnow

% bivariate relationships
axes(axh(5));
imagesc(mbyc, [-1 1]); %title('Component'); ylabel('Movement');
set(gca,'YLim', [1 n], 'XLim', [.5 m+.5], 'YDir', 'Reverse');
axis off
axes(axh(6));
imagesc(mbyt, [-1 1]); %title('Task'); 
set(gca,'YLim', [1 n], 'XLim', [.5 p+.5], 'YDir', 'Reverse');
axis off
axes(axh(8));
imagesc(tbyc, [-1 1]); xlabel('Component'); %ylabel('Task');
set(gca,'YLim', [1 p], 'XLim', [.5 m+.5], 'YDir', 'Reverse');
set(gca,'YColor',[1 1 1])

axes(axh(9));
axis off


% multivariate relationships
h3 = axes('position',[0.4108    0.7093-.08    0.2134    0.01]);
rsq_mvmt = rsquare_calc(rpz, V);
imagesc(rsq_mvmt, [-1 1]);
title('R^2 by Mvmt', 'FontSize', 14);
axis off

h4 = axes('position',[0.4108    0.4096-.08    0.2134    0.01]);
rsq_task = rsquare_calc(X, V);
imagesc(rsq_task, [-1 1]);
title('R^2 by Task', 'FontSize', 14);
axis off

h4 = axes('position',[0.6916    0.7093-.08    0.2134    0.01]);
rsq_mvmt = rsquare_calc(rpz, X);
imagesc(rsq_mvmt, [-1 1]);
title('R^2 by Mvmt', 'FontSize', 14);
axis off

% axis colorbar
h = axes('Position', [0.11    0.7093 .2 .01]);
imagesc([0], [-1 1]); colorbar('horiz');
title('Correlation', 'FontSize', 14);
axis off

% Var Inflation
h2 = axes('Position', [.69 .27 .2134 .04], 'FontSize', 14, 'YColor', [1 1 1]);
imagesc(vifs, [0 10]);
%colorbar('horiz')
title('Var. Inflation (VIF)');
axis off

% Var Inflation
h2 = axes('Position', [.69 .19 .2134 .04], 'FontSize', 14, 'YColor', [1 1 1]);
imagesc(vifsm, [0 10]);
colorbar('horiz')
title('VIF w/Mvmt');

%cmap = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
cmap = colormap_tor([0 0 1], [1 1 0], [.5 .2 1], [1 1 1], [1 0 0]);
colormap(cmap)


print('-dpsc2', '-append', 'qc_report');
scn_export_papersetup(800); saveas(gcf,['qc_images' filesep 'task_mvmt_component_summary'],'png');

end

