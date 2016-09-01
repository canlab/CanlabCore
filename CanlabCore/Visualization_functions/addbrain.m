function p = addbrain(varargin)
% :Usage:
% ::
%
%    handle = addbrain([method],enter 2nd arg to suppress lighting changes)
%
% quick function to add transparent brain surface to figure
% ::
%
%    han = addbrain;   % lateral surface
%
%    han = addbrain('brainstem');
%
% NOTE: this version uses structures in SPM2 space (Colin atlas)
% Available keywords:
%
% :CORTICAL SURFACES:
%
%   **'transparent_surface':**
%        the default.  2 mm res SPM2 brain surface
%
%   **'hires':**
%        a high-resolution surface (from Caret segmentation)
%
%   **'hires left':**
%        hi-resolution left medial with cerebellum (Caret seg)
%
%   **'hires right':**
%        same, right hem
%
%   **'surface left':**
%        hi-resolution left medial with cerebellum
%        slightly expanded, from Glasser et al. 2016 Nature 
%        based on the Human Connectome Project
%
%   **'surface right':**
%        same, right hem
%
%   **'left':**
%        2 mm resolution left hem, no cerebellum
%
%   **'right':**
%
%   **'vmpfc':**
%
% :CUTAWAY SURFACES:
%
%   **'brainbottom':**
%
% :COMPOSITES:
%
%   **'limbic':**
%        A collection of subcortical nuclei with left surface
%
%   **'foursurfaces':**
%        Lateral and medial views, with brainstem
%
%   **'BG':**
%        Basal ganglia
%
% :SUBCORTICAL SURFACES:
%   - 'brainstem'
%   - 'suit brainstem'
%   - 'amygdala'
%   - 'thalamus'
%   - 'hippocampus'
%   - 'midbrain'
%   - 'caudate'
%   - 'globus pallidus'
%   - 'putamen'
%   - 'nucleus accumbens'
%   - 'hypothalamus'
%   - 'cerebellum'
%   - {'md','mediodorsal'}
%   - {'cm','centromedian'}
%   - 'pbn'
%   - 'rvm'
%   - 'nts'
%   - 'lc'
%   - {'sn', 'substantia nigra'}
%   - {'stn', 'subthalamic nucleus'}
%   - {'rn', 'red nucleus'}
%   - {'olive', 'inferior olive'}
%   - {'nrm', 'raphe magnus'}
%
% :SPECIAL COMMANDS:
% 
% han = addbrain('colorchange',my_rgb_color,han);
%
% Change gray background to some other color, excluding blobs already rendered
%  - han: Input handles with patch object
%  - my_rgb_color: [x x x] color triplet
%  - Only works for changing from gray background right now. 
%
% han = addbrain('eraseblobs',han);
%
% Set all rendered blob colors back to gray; useful for re-rendering on existing surfaces.
%  - han: Input handles with patch object
%  - Only works for changing to gray background right now. 
%
% See also: cluster_surf, img2surf.m, surface() methods for objects, cluster_cutaways

p = [];
meth = 'transparent_surface';

docolor = 1;
if length(varargin) > 0
    meth = varargin{1};
end

if length(varargin) > 1
    color = varargin{2};
end

switch meth

    case 'colorchange'
        myp = varargin{3};
        if isempty(myp), myp = findobj(gcf,'Type','Patch'); end
        background_colorchange(myp,color);
        docolor = 0;
        
    case {'eraseblobs', 'set_all_to_gray'}
        myp = varargin{2};
        for i = 1:length(myp)
            %set(hh, 'FaceColor', [.5 .5 .5]);
            len = size(myp(i).FaceVertexCData, 1);
            set(myp(i), 'FaceVertexCData', repmat([.5 .5 .5], len, 1));
        end
        p = myp;
        return
        
    case 'left'
        pname = 'surf_spm2_left.mat';  % moderate res, no cerebellum

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);

        view(90,0); axis off; axis image; lightRestoreSingle(gca); material dull;

    case 'hires left'
        pname = 'surf_spm2_brain_left.mat'; % high res, with cblm.  caret segmentation

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);

        view(90,0);
        axis off; 
        axis image; 
        lightRestoreSingle(gca); 
        material dull;
        
    case 'surface left'
       pname = 'L.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat'; % from Glasser_et_al_2016_HCP
        p1 = add_surface(pname);
        set(p1,'FaceColor',[.5 .5 .5]);
         p2 = add_surface('suit_surface_brainstem_cerebellum.mat');
        set(p2,'FaceColor',[.5 .5 .5]);
        view(90,0);
        axis off;
        axis image;
        %         lightRestoreSingle(gca);
        material dull;
        p=[p1 p2];
        
    case 'right'

        pname = 'surf_spm2_right.mat'; %'surf_single_subj_grayR.mat';

        p = add_surface(pname);

        set(p,'FaceColor',[.5 .5 .5]);
        view(270,0); axis off; axis image; lightRestoreSingle(gca); material dull;

    case 'hires right'
        pname = 'surf_spm2_brain_right.mat'; % high res, with cblm.  caret segmentation

        p = add_surface(pname);

        set(p,'FaceColor',[.5 .5 .5]);
        view(270,0); 
        axis off; 
        axis image; 
        lightRestoreSingle(gca);
        material dull;
        
    case 'surface right'
        pname = 'R.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat'; % from Glasser_et_al_2016_HCP
        p1 = add_surface(pname);
        set(p1,'FaceColor',[.5 .5 .5]);
         p2 = add_surface('suit_surface_brainstem_cerebellum.mat');
        set(p2,'FaceColor',[.5 .5 .5]);
        view(90,0);
        axis off;
        axis image;
        %         lightRestoreSingle(gca);
        material dull;
        p=[p1 p2];
        
    case 'transparent_surface'

        %spm99 pname = 'surf_single_subj_T1_gray.mat';  %'surf_single_subj_gw_sparse.mat'; %
        pname = 'surf_spm2_brain.mat';  % medium res, caret segmentation
%         pname = 'LR.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat'; % from Glasser_et_al_2016_HCP

        p = add_surface(pname);

        set(p, 'FaceAlpha', .7)
        lighting gouraud
        lightRestoreSingle
        axis vis3d image tight
        
    case 'hires'

        pname = 'surf_spm2_brain_1mm.mat';  % hi res, caret segmentation

        p = add_surface(pname);

        set(p, 'FaceAlpha', .7)
        lighting gouraud
        lightRestoreSingle;
        axis vis3d image tight
        
        
    case 'brainstem'

        pname = 'surf_spm2_brainstem.mat';
%         pname = 'suit_surface_brainstem_cerebellum.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .65 .4]);
        
    case 'suit brainstem'

        pname = 'suit_surface_brainstem_cerebellum.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);
         axis off;
        axis image;
        material dull;
        
    case 'brainbottom'
        [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -20],'filename','scalped_single_subj_T1');
        set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
        h = findobj('Type','Light'); delete(h); [az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
        set(p,'FaceAlpha',1)

    case 'amygdala'
        % Carmack
        %         P = which('Tal_Amy.img');
        %         [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
        %         if findstr(P,'Tal_Amy.img'), str(49:58) = '[.4 .4 .6]'; , end

        % ICBM
        %             load amy_clusters
        %             p = imageCluster('cluster',amy,'color',[0 0 .5],'alpha',.5);

        pname = 'surf_spm2_amy.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[0 0 .5]);

    case 'thalamus'
        %         P = which('carmack_thal_bstem.mat'); load(P)
        %         p = imageCluster('cluster',thal,'color',[0 .8 .3],'alpha',.5);
        pname = 'surf_spm2_thal.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.9 .65 .5]);

    case {'hippocampus', 'hipp'}
        %Carmack
        %                      P = which('Tal_Hip.img');
        %                      [p,outP,FV, cl, myLight] = mask2surface(P,0,[.5 .5 0]);
        %                      if findstr(P,'Tal_Hip.img'), str(49:58) = '[.5 .6 .6]';, end

        % ICBM
        %load hipp_clusters hipp
        %p = imageCluster('cluster',hipp,'color',[.5 .6 .6],'alpha',.5);

        pname = 'surf_spm2_hipp.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.7 .4 .4]);

    case 'midbrain'
        P = which('carmack_thal_bstem.mat');
        load(P)
        p = imageCluster('cluster',midbrain,'color',[.7 .3 0],'alpha',.5);

    case 'caudate'
        %P = which('Tal_Cau.img'); %which('ICBM_caudate.img');   %('Tal_Cau.img'); %
        %[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
        %if findstr(P,'Tal_Cau.img'), str(49:58) = '[.5 .6 .6]'; delete(p(1)); p = p(2);,eval(str),end
        %         P = which('carmack_more_clusters.mat'); load(P)
        %         p = imageCluster('cluster',cau,'color',[.5 .6 .6],'alpha',.5);

        pname = 'surf_spm2_caudate.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .6 .6]);

    case {'globus pallidus', 'gp'}
%         P = which('Tal_Glo.img');
%         [p,outP,FV, cl, myLight] = mask2surface(P,0,[.5 .6 .5]);
%         if findstr(P,'Tal_Glo.img'), str(49:58) = '[.5 .6 .5]'; , end
        
        load(which('Keuken_2014_7T_regions.mat'), 'GPe', 'GPi')
        p = imageCluster('cluster', GPe, 'color',[.5 .6 .5],'alpha', .5);
        set(p(end), 'Tag', 'GPe');
        p = [p imageCluster('cluster', GPi, 'color',[.5 .75 .5],'alpha', .5)];
        set(p(end), 'Tag', 'GPi');
        
    case {'putamen', 'put'}
        %P = which('Tal_Put.img');
        %[p,outP,FV, cl, myLight] = mask2surface(P,0,[.9 .4 0]);
        %if findstr(P,'Tal_Put.img'), str(49:58) = '[.5 .5 .6]'; , end
        %P = which('carmack_more_clusters.mat'); load(P)]
        
%         fbase = 'LBPA40_spm5_label_clusters.mat';
%         fname = which(fbase); if isempty(fname), error(['Looking for ' fbase]); end
%         load(fname)
%         put = cat(2, cl{51:52});
%         p = imageCluster('cluster',put,'color',[.5 .5 .6],'alpha',.5);
%         
        pname = 'spm_surf_putamen_luke.mat';
         p = add_surface(pname);
        set(p,'FaceColor',[.3 .7 .5]);
                   
    case {'BG', 'bg', 'basal ganglia'}
        p = addbrain('gp');
        
        p = [p  addbrain('put')];
        set(p(end), 'Tag', 'putamen');
        
        p = [p  addbrain('caudate')];
        set(p(end), 'Tag', 'caudate');
         
         
    case {'nucleus accumbens','nacc','nac'}
        %             P = which('NucAccumb_clusters.mat');
        %             load(P)
        %             p(1) = imageCluster('cluster',cl(3),'color',[0 .5 0],'alpha',.5);
        %             p(2) = imageCluster('cluster',cl(4),'color',[0 .5 0],'alpha',.5); ,

        pname = 'surf_spm2_nac.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[0 .5 0]);

    case {'hypothalamus','hy','hythal'}
        %P = which('Tal_Put.img');
        %[p,outP,FV, cl, myLight] = mask2surface(P,0,[.9 .4 0]);
        %if findstr(P,'Tal_Put.img'), str(49:58) = '[.5 .5 .6]'; , end

        %         P = which('carmack_more_clusters.mat'); load(P)
        %         p = imageCluster('cluster',hy,'color',[1 1 0],'alpha',.5);

        % Carmack
        %         load hy_clusters hy
        %         p = imageCluster('cluster',hy,'color',[1 1 0],'alpha',.5);

        pname = 'surf_spm2_hythal.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.9 .85 .5]);

    case {'cerebellum','cblm'}
        pname = 'surf_spm2_cblm.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[.8 .65 .8]);

    case 'limbic'
        p = [];
        myp = addbrain('amygdala');p = [p myp];
        myp = addbrain('hypothalamus');p = [p myp];
        myp = addbrain('hippocampus');p = [p myp];
        myp = addbrain('thalamus');p = [p myp];
        myp = addbrain('nacc');p = [p myp];
        myp = addbrain('hires left');p = [p myp];
         myp = addbrain('BG');p = [p myp];
%          myp = addbrain('caudate');p = [p myp];
%         myp = addbrain('gp'); p = [p myp];
%          myp = addbrain('putamen'); p = [p myp];
         myp = addbrain('brainstem'); p = [p myp];
        set(p,'FaceAlpha',1);

        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)


    case 'pag'
        load pag_cl
        p = imageCluster('cluster',pag,'color',[1 0 0],'alpha',1);

        %P = which('ROI_pag.img');
        %P = which('spm5_pag.img');
        %[p,outP,FV, cl, myLight] = mask2surface(P,0,[1 0 0]);

    case {'md','mediodorsal'}
        load thal_brainstem_approx_working MD
        p = imageCluster('cluster',MD,'color',[1 0 0],'alpha',1);

    case {'cm','centromedian'}
        load thal_brainstem_approx_working CM
        p = imageCluster('cluster',CM,'color',[1 .7 0],'alpha',1);

    case 'pbn'
        load pbn_cl
        p = imageCluster('cluster',pbn,'color',[1 .5 0],'alpha',1);

    case 'rvm'
        load rvm_cl
        p = imageCluster('cluster',rvm,'color',[1 .2 .1],'alpha',1);

    case 'nts'
        load nts_cl
        p = imageCluster('cluster',nts,'color',[0 0 1],'alpha',1);
    
    case {'lc'} %******
        P = which('ROI_LC.img');
        [p,outP,FV, cl, myLight] = mask2surface(P,0,[1 .5 0]);

    case {'sn', 'substantia nigra'}
%         P = which('ROI_SN.img');
%         [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);

        load(which('Keuken_2014_7T_regions.mat'), 'SN')
        p = imageCluster('cluster', SN, 'color',[0 0 .5],'alpha', .5);
        
    case {'stn', 'subthalamic nucleus'} % ****
%         P = which('ROI_STN.img');
%         [p,outP,FV, cl, myLight] = mask2surface(P,0,[1 0 .5]);

        load(which('Keuken_2014_7T_regions.mat'), 'STN')
        p = imageCluster('cluster', STN, 'color',[1 0 .5],'alpha', .5);
        
    case {'rn', 'red nucleus'}
%         P = which('ROI_red_nucleus.img');
%         [p,outP,FV, cl, myLight] = mask2surface(P,0,[1 0 0]);

        load(which('Keuken_2014_7T_regions.mat'), 'RN')
        p = imageCluster('cluster', RN, 'color',[1 0 0],'alpha', .5);
        
    case {'olive', 'inferior olive'}
        P = which('ROI_inf_olive.img');
        [p,outP,FV, cl, myLight] = mask2surface(P,0,[.5 1 .5]);

    case {'nrm', 'raphe magnus'}
        P = which('ROI_raphe_magnus.img');
        [p,outP,FV, cl, myLight] = mask2surface(P,0,[.3 0 1]);

    case {'vmpfc', 'vmPFC'}
        
        P = which('VMPFC_display_mask.img');
         [p,outP,FV, cl, myLight] = mask2surface(P, 0, [.7 .3 0]);

    case 'foursurfaces'
        P = run_foursurfaces;
        
    case 'flat left'
        
        pname = which('L.flat.32k.mat'); % from Glasser_et_al_2016_HCP
        
        p = add_surface(pname);
        
        view(0,90); 
        axis off; 
        axis image; 
%         material dull;
        
    case 'flat right'
        
        pname = which('R.flat.32k.mat'); % from Glasser_et_al_2016_HCP
        
        p = add_surface(pname);
        
        view(0,90); 
        axis off; 
        axis image; 
%         material dull;
  
    otherwise
        error('Unknown method.');
        
end  % method


switch meth
    
    case {'limbic', 'BG', 'globus pallidus', 'gp', 'bg', 'basal ganglia'}
        % do nothing; multi-region
        
    otherwise
        if all(ishandle(p))
            set(p, 'Tag', meth);
        end
end

if docolor && exist('color','var')
    if isa(color,'double')
    set(p,'FaceColor',color);
    elseif isa(color,'patch')
        
    end
end


% suppress lighting if 2nd arg, otherwise, do it.

if length(varargin) > 2
    lightRestoreSingle(gca); 
    camlight right;
    %set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView
end

lighting gouraud
axis vis3d image tight
material dull
drawnow

%view(135,30)


end % main function




function p = add_surface(pname)
Ps = which(pname); %'c:\tor_scripts\3DheadUtility\surf_single_subj_T1_gray.mat';
if isempty(Ps), disp(['I need the file: ' pname]); return; end

%Ps = which('surf_single_subj_grayR.mat');
%Ps = which('surf_brain_render_T1_preCarmack.mat');
load(Ps)

if exist('cdata','var')
%     p = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',cdata, ...
%         'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200,'FaceColor','interp');
%     colormap(gray);
 p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
    set(p,'FaceAlpha',.3)
    
else
    p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
    set(p,'FaceAlpha',.3)
    
end


end


function background_colorchange(myp,mycolor)


for i = 1:length(myp)
    p = myp(i);


    % change non-active to input color
    vdat = get(p,'FaceVertexCData');
    wh = find(all(vdat == .5,2));
    vdat(wh,:) = repmat(mycolor,length(wh),1);
    set(p,'FaceVertexCData',vdat);

end

end


% Four surface subfunction
% ---------------------------------------------------------------------
function all_surf_handles = run_foursurfaces

all_surf_handles = [];
f1 = gcf;

nrows = 2;
ncols = 2;

% Right lateral
% ------------------------------------------------------------------------
figure(f1);
subplot(nrows, ncols , 1);
surfh = addbrain('hires right');
set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
view(90, 0)
lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

surfh2 = addbrain('brainstem');
surfh2 = [surfh2 addbrain('thalamus')];
set(surfh2, 'FaceColor', [.5 .5 .5], 'FaceAlpha', .8);
surfh = [surfh surfh2];

all_surf_handles = [all_surf_handles surfh];

% Right medial
% ------------------------------------------------------------------------

figure(f1);
axh = subplot(nrows, ncols , 4);
copyobj(surfh, axh);
view(270, 0);

lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

% Left lateral
% ------------------------------------------------------------------------

figure(f1);
subplot(nrows, ncols , 2);
surfh = addbrain('hires left');
set(surfh, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 1);
view(270, 0)
lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

surfh2 = addbrain('brainstem');
surfh2 = [surfh2 addbrain('thalamus')];
set(surfh2, 'FaceColor', [.5 .5 .5], 'FaceAlpha', .8);
surfh = [surfh surfh2];

all_surf_handles = [all_surf_handles surfh];

% Left medial
% ------------------------------------------------------------------------

figure(f1);
axh = subplot(nrows, ncols , 3);
copyobj(surfh, axh);
view(90, 0);

lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

end

