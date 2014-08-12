function OUT = xval_regression_bootstrapweightmap_brainplots(STATS, varargin)
% OUT = xval_regression_bootstrapweightmap_brainplots(STATS, varargin)
%
% Creates plots for xval regression bootstrap STATS structure
% See xval_regression_multisubject_bootstrapweightmap.m
%
% Options:
% =================================================
% {'newmontage', 'montage'}, donewmontage = 1;
% {'oldmontage'}, dooldmontage = 1;
% {'surface'}, dosurface = 1;
% {'diary'}, dodiary = 1;
% {'save'}, dosave = 1;
% 'reversecolors', reversecolors = 1;
%
% OUT has the thresholds and positive/negative clusters for each threshold
%
% Examples:
% =================================================
% This saves all output except old-style montages in the current directory:
% xval_regression_bootstrapweightmap_brainplots(STATS,'montage', 'surface', 'diary', 'save');

                
donewmontage = 0;
dooldmontage = 0;
dosurface = 0;
dodiary = 0;
dosave = 0;
reversecolors = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'newmontage', 'montage'}, donewmontage = 1;
                
            case {'oldmontage'}, dooldmontage = 1;
            case {'surface'}, dosurface = 1;
            case {'diary'}, dodiary = 1;
            case {'save'}, dosave = 1;
                
            case 'reversecolors', reversecolors = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


volInfo = STATS.VOXWEIGHTS.volInfo;
w = STATS.VOXWEIGHTS.vox_weights{1};
p = STATS.VOXWEIGHTS.pval_2tail; 

% pvals will be min for out-of-mask voxels
p(w' == 0) = 1;

% Print summary table
% ----------------------------------------------

if dodiary
    
    diary bootstrap_weightmap_output.txt
    
end

% Thresholds
% ----------------------------------------------

thresh = [Inf .001 .002 .005 .01 .05];
sigvox = sum(STATS.VOXWEIGHTS.sigfdr);
threshnames = {'FDR' '001' '002' '005' '01' '05'};

thresh(1) = STATS.VOXWEIGHTS.pthr_fdr;

% Sig voxels and clusters
% ----------------------------------------------

for ii = 2:length(thresh)
    sigvox(ii) = sum(p < thresh(ii) & w' ~= 0);
end

fprintf('In-Mask Space is %3.0f voxels\n', volInfo.n_inmask)
fprintf('Valid non-empty voxels (for FDR calculation): %3.0f voxels\n\n', sum(w' ~= 0))

fprintf('FDR threshold = %3.6f, sig vox = %3.0f\n', STATS.VOXWEIGHTS.pthr_fdr, sum(STATS.VOXWEIGHTS.sigfdr))
for ii = 2:length(thresh)
    fprintf('p < %3.4f, sig vox = %3.0f\n', thresh(ii), sigvox(ii))
end

if dodiary
    
    diary off
    
end

OUT = struct('thresh', thresh, 'sigvox', sigvox, 'threshnames', threshnames, 'clpos', [], 'clneg', []);

% Set Colors
% ----------------------------------------------

if reversecolors, revstr = '_revcolors'; else revstr = ''; end

if reversecolors
    negcm = colormap_tor([.8 0 0], [1 1 0], [1 1 0]);
    poscm = colormap_tor([.4 .8 .6], [0 0 1], [0 0 1]);
else
    poscm = colormap_tor([.8 0 0], [1 1 0], [1 1 0]);
    negcm = colormap_tor([.4 .8 .6], [0 0 1], [0 0 1]);
end

% Save weight image
% ----------------------------------------------
w = STATS.VOXWEIGHTS.vox_weights{1};

if dosave
    iimg_reconstruct_vols(w, volInfo, 'outname', 'voxelweights.img');
    
    iimg_reconstruct_vols(p', volInfo, 'outname', 'voxelweight_pvalues.img');
    
end

% -------------------------------------------------------------
% Save weight images and display montages for each threshold
% -------------------------------------------------------------

% NEW MONTAGE SETUP
% ----------------------------------------------
if donewmontage
    obj = fmridisplay;            % Create object with canonical underlay
    obj = montage(obj, 'onerow');  % Show axial montage of underlay
    enlarge_axes(gcf, .95);
end

for i = 1:length(thresh)
    
    if ~sigvox(i), continue, end  % skip if no vox
    
    outbase = sprintf('voxelweights_p_%s%s', threshnames{i}, revstr);
    
    montagebase = sprintf('thresh_p_%s_%3.0fvox%s', threshnames{i}, sigvox(i), revstr);
    mb = montagebase; mb(mb == '_') = ' ';
    
    disp(outbase)
    
    % Save weight image
    % ----------------------------------------------
    outname = [outbase '.img'];
    wthr = w;
    wthr(p > thresh(i)) = 0;
    
    if dosave
        iimg_reconstruct_vols(wthr, volInfo, 'outname', outname);
    end
    
    % ***Note: this can fail if no variability in Z values, i.e. if
    % few booot samples are used. should fix...***
    
    % Clusters
    zthr = sign(w) .* norminv(1 - p');
    zthr(p > thresh(i)) = 0;
    %             wpos = wthr; wpos(wpos < 0) = 0;
    %             wneg = wthr; wneg(wneg > 0) = 0;
    zpos = zthr; zpos(zpos < 0) = 0;
    zneg = zthr; zneg(zneg > 0) = 0;
    clpos = iimg_indx2clusters(zpos, volInfo);
    clneg = iimg_indx2clusters(zneg, volInfo);
    clpos(cat(1, clpos.numVox) < 3) = [];
    clneg(cat(1, clneg.numVox) < 3) = [];
        
    OUT(i).clpos = clpos;
    OUT(i).clneg = clneg;
    
    % NEW MONTAGE blobs
    % -------------------------------
    if donewmontage
        
        obj = removeblobs(obj);

        obj = addblobs(obj, clpos, 'maxcolor', [1 1 0], 'mincolor', [1 .3 0]);
        obj = addblobs(obj, clneg, 'maxcolor', [.5 0 1], 'mincolor', [0 0 1]);
        
        axes(obj.montage{1}.axis_handles(2))
        title([mb ' k=3 contig'], 'FontSize', 18);
        
        obj = legend(obj);
        
        if dosave
            scn_export_papersetup(500); 
            saveas(gcf, [montagebase '_new_axial.png']);
            saveas(gcf, [montagebase '_new_axial.fig']);
        end
        
    end
    
    if dosave
        diary bootstrap_weightmap_output.txt
    end
    
    fprintf('\n===============================\n')
    
    disp([mb 'k=3 contig'])
    fprintf('===============================\n')
    fprintf('POSITIVE EFFECTS\n_______________________________\n')
    
    cluster_table(clpos, 0 , 0);
    
    fprintf('\nNEGATIVE EFFECTS\n_______________________________\n')
    
    cluster_table(clneg, 0 , 0);
    
    if dooldmontage
        
        % Orthviews
        % ----------------------------------------------
        cl = iimg_indx2clusters(wthr, volInfo);
        cluster_orthviews(cl)
        spm_orthviews_hotcool_colormap(w, [prctile(w, 50)]);
        cm = get(gcf, 'Colormap');
        
        % adjust color map
        %isgray = all(cm - repmat(mean(cm, 2), 1, 3)  < .01, 2)
        cm2 = [cm(1:64, :); negcm(end:-2:1, :); poscm(1:2:end, :)];
        colormap(cm2)
        
        % Montages from orthviews
        % ----------------------------------------------
        overlay = which('SPM8_colin27T1_seg.img');
        
        slices_fig_h = cluster_orthviews_montage(6, 'axial', overlay, 0, 'onerow');
        %colormap(cm); h = findobj(gcf, 'Type', 'Text'); delete(h)
        if dosave
            scn_export_papersetup(500); saveas(gcf, [montagebase '_axial.png']);
        end
        
        slices_fig_h = cluster_orthviews_montage(6, 'sagittal', overlay, 0, 'onerow');
        %colormap(cm); h = findobj(gcf, 'Type', 'Text'); delete(h)
        if dosave
            scn_export_papersetup(500); saveas(gcf, [montagebase '_sagittal.png']);
        end
        
        slices_fig_h = cluster_orthviews_montage(6, 'coronal', overlay, 0, 'onerow');
        %colormap(cm); h = findobj(gcf, 'Type', 'Text'); delete(h)
        if dosave
            scn_export_papersetup(500); saveas(gcf, [montagebase '_coronal.png']);
        end
        
    end
    
end


% -----------------------------------------------
% Figure - surface, varying threshold
% -----------------------------------------------

if dosurface
    
    for i = 1:length(thresh)
        
        if ~sigvox(i), continue, end  % skip if no vox
        
        outname = sprintf('thresh_p_%s_%3.0fvox%s', threshnames{i}, sigvox(i), revstr);
        
        disp(outname)
        
        %outbase = 'voxelweights_p002';
        %outname = [outbase '.img'];
        
        wthr = w; wthr(p > thresh(i)+eps) = 0;
        cl = iimg_indx2clusters(wthr, volInfo);
        
        create_figure('Brain_Surface', 2, 2);
        
        
        han1 = addbrain('hires');
        set(han1, 'FaceAlpha', 1);
        cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han1);
        
        %replace with 'hires' and re-run to do whole cortical surface
        create_figure('Brain_Surface', 2, 2, 1);
        subplot(2, 2, 2);
        han = addbrain('hires right'); set(han, 'FaceAlpha', 1);
        
        drawnow
        cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
        set(han, 'FaceAlpha', 1);
        axis image; lightRestoreSingle;
        lighting gouraud
        
        create_figure('Brain_Surface', 2, 2, 1);
        subplot(2, 2, 3);
        han = addbrain('hires left'); set(han, 'FaceAlpha', 1);
        axis image; lightRestoreSingle;
        lighting gouraud
        drawnow
        cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
        
        create_figure('Brain_Surface', 2, 2, 1);
        subplot(2, 2, 4);
        han = addbrain('limbic'); set(han, 'FaceAlpha', 1);
        axis image; lightRestoreSingle;
        lighting gouraud
        drawnow
        cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
        
        set(han(end), 'FaceAlpha', .2)
        han2 = addbrain('brainstem');
        cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, han);
        view(135, 15)
        
        create_figure('Brain_Surface', 2, 2, 1);
        subplot(2, 2, 1)
        han = findobj(gca,'Type','patch'); set(han, 'FaceAlpha',1)
        axis image
        
        if dosave
            scn_export_papersetup(600); saveas(gcf, [outname '_Surface1'], 'png');
            subplot(2, 2, 1);
            view(135, 20); lightRestoreSingle;
            saveas(gcf, [outname '_Surface2'], 'png');
            view(225, 20); lightRestoreSingle;
            saveas(gcf, [outname '_Surface3'], 'png');
            view(90, 3); lightRestoreSingle;
            subplot(2, 2, 4);
            view(135, 0); lightRestoreSingle;
            saveas(gcf, [outname '_Surface4'], 'png');
            subplot(2, 2, 1);
            view(270, 3); lightRestoreSingle;
            subplot(2, 2, 4);
            view(235, 0); lightRestoreSingle;
            saveas(gcf, [outname '_Surface5'], 'png');
            subplot(2, 2, 1);
            view(180, -90); lightRestoreSingle;
            saveas(gcf, [outname '_Surface6'], 'png');
        end
        
    end % threshold loop
    
end % surface


end