function p = addbrain(varargin)
% handle = addbrain([method],enter 2nd arg to suppress lighting changes)
% quick function to add transparent brain surface to figure
%
% % han = addbrain;   % lateral surface
% han = addbrain('brainstem');
%
% NOTE: this version uses structures in SPM2 space (Colin atlas)
% Available keywords:
% 'transparent_surface' : the default.  2 mm res SPM2 brain surface
% 'hires'       : a high-resolution surface (from Caret segmentation)
% 'hires left'  : hi-resolution left medial with cerebellum (Caret seg)
% 'hires right' : same, right hem
% 'left'        :  2 mm resolution left hem, no cerebellum
% 'right'
% 'brainstem'
% 'brainbottom'
% 'amygdala'
% 'thalamus'
% 'hippocampus'
% midbrain'
%     'caudate'
%     'globus pallidus'
%     'putamen'
%     'nucleus accumbens'
%     'hypothalamus'
%     'cerebellum'
% case {'md','mediodorsal'}
%     case {'cm','centromedian'}
%     case 'pbn'
%     case 'rvm'
%     case 'nts'
%     case {'lc'}
%     case {'sn', 'substantia nigra'}
%     case {'stn', 'subthalamic nucleus'}
%     case {'rn', 'red nucleus'}
%     case {'olive', 'inferior olive'}
%     case {'nrm', 'raphe magnus'}
%
% han = addbrain('colorchange',my_rgb_color,han);
% Works on patch object in input handle han
% or, if han is empty, on all patch objects in the current figure.
% Only works for changing from gray background right now.  (Future update to handle any
% background.)

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

    case 'left'
        pname = 'surf_spm2_left.mat';  % moderate res, no cerebellum

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);

        view(90,0); axis off; axis image; lightRestoreSingle(gca); material dull;

    case 'hires left'
        pname = 'surf_spm2_brain_left.mat'; % high res, with cblm.  caret segmentation

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);

        view(90,0); axis off; axis image; lightRestoreSingle(gca); material dull;

    case 'right'

        pname = 'surf_spm2_right.mat'; %'surf_single_subj_grayR.mat';

        p = add_surface(pname);

        set(p,'FaceColor',[.5 .5 .5]);
        view(270,0); axis off; axis image; lightRestoreSingle(gca); material dull;

    case 'hires right'
        pname = 'surf_spm2_brain_right.mat'; % high res, with cblm.  caret segmentation

        p = add_surface(pname);

        set(p,'FaceColor',[.5 .5 .5]);
        view(270,0); axis off; axis image; lightRestoreSingle(gca); material dull;

    case 'transparent_surface'

        %spm99 pname = 'surf_single_subj_T1_gray.mat';  %'surf_single_subj_gw_sparse.mat'; %
        pname = 'surf_spm2_brain.mat';  % medium res, caret segmentation

        p = add_surface(pname);

    case 'hires'

        pname = 'surf_spm2_brain_1mm.mat';  % hi res, caret segmentation

        p = add_surface(pname);

    case 'brainstem'

        pname = 'surf_spm2_brainstem.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.5 .65 .4]);

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
        p = [p imageCluster('cluster', GPi, 'color',[.5 .75 .5],'alpha', .5)];
        
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
        set(p,'FaceColor',[.5 .9 .6]);
        
    case {'BG', 'bg', 'basal ganglia'}
        p = addbrain('gp');
        p = [p  addbrain('put')];
        p = [p  addbrain('caudate')];

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
        myp = addbrain('caudate');p = [p myp];
        myp = addbrain('left');p = [p myp];
        myp = addbrain('gp'); p = [p myp];
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


    otherwise
        error('Unknown method.');
end


switch meth
    
    case 'limbic'
        % do nothing; multi-region
    otherwise
        if all(ishandle(p))
            set(p, 'Tag', meth);
        end
end

if docolor && exist('color','var')
    set(p,'FaceColor',color);
end


% suppress lighting if 2nd arg, otherwise, do it.

if length(varargin) > 2
    lighting gouraud;
    axis image; %myLight = camlight(0,0);set(myLight,'Tag','myLight');
    lightRestoreSingle(gca); camlight right
    %set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView
end


drawnow
axis vis3d
%view(135,30)


return




function p = add_surface(pname)
Ps = which(pname); %'c:\tor_scripts\3DheadUtility\surf_single_subj_T1_gray.mat';
if isempty(Ps), disp(['I need the file: ' pname]); return; end

%Ps = which('surf_single_subj_grayR.mat');
%Ps = which('surf_brain_render_T1_preCarmack.mat');
load(Ps)
p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
    'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);

set(p,'FaceAlpha',.3)
return


function background_colorchange(myp,mycolor)


for i = 1:length(myp)
    p = myp(i);


    % change non-active to input color
    vdat = get(p,'FaceVertexCData');
    wh = find(all(vdat == .5,2));
    vdat(wh,:) = repmat(mycolor,length(wh),1);
    set(p,'FaceVertexCData',vdat);

end

return
