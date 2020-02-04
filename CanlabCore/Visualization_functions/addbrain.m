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
% NOTE: this version uses structures in MNI space
% The exact space depends on the method used.
%
% % Regions
% % -----------------------------------------------------------------------
% 'vmpfc' 'nacc' 'BST' 'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
% 'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
% 'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
% 'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'sc' 'ic' 'drn' 'mrn' ...
% 'thalamus' 'thal' 'LGN' 'lgn' 'MGN' 'mgn' 'VPthal', 'VPLthal', 'VPL', 'intralaminar_thal', ...
% 'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'pbn' 'ncf' 'vep' 'PBP'
% 'amygdala' 'amygdala hires' 'hippocampus', 'hipp' 'hippocampus hires' 
% 
% % Cortical Surfaces
% % -----------------------------------------------------------------------
% 'left' 'hires left' 'surface left' 'hires surface left' ...
% 'right' 'hires right' 'surface right' 'hires surface right' ...
% 'transparent_surface' 'foursurfaces' 'flat left'  'flat right' ...
% 
% % Macro subcortical surfaces
% % -----------------------------------------------------------------------
% 'CIT168' 'cerebellum','cblm' 'brainstem' 'suit brainstem' 
% 
% % Cutaways
% % -----------------------------------------------------------------------
% 'brainbottom' 'cutaway', 'left_cutaway' 'right_cutaway' ...
% 'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5' ...
% 
% % Groups
% % -----------------------------------------------------------------------
% 'bg', 'basal ganglia' 'midbrain_group' 'limbic' 'limbic hires' 'brainstem_group' 'thalamus_group' 
%
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
%       Bottom of brain and head
%
%   **'cutaway':**
%       A canonical surface cutaway 
%       Uses: canlab_canonical_brain_surface_cutaways
%       pre-2020 used surface_cutaway.m
%
%   **'left_cutaway'**
%   **'right_cutaway' **
%   **'left_insula_slab'**
%   **'right_insula_slab'**
%   **'accumbens_slab'**
%
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
%   **'midbrain_group'**
%        Midbrain structures
%
%   **'brainstem_group'**
%        Midbrain and pons/medulla structures
%
%   **'thalamus_group'**
%        Midbrain and pons/medulla structures
%
% :SUBCORTICAL SURFACES: see canlab_load_ROI.m for complete list
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
% Examples:
% figure; addbrain('midbrain_group');
% addbrain('lc'); addbrain('rvm'); addbrain('VPL'); addbrain('thalamus');
% addbrain('bg');
% addbrain('hires left');
% view(135, 10); lightRestoreSingle;
%
% See also: canlab_load_ROI, cluster_surf, img2surf.m, surface() methods for objects, cluster_cutaways

%%% Programmer's notes:
%
% 2017/03/03 added Jon Brooks' RVM and moved call to old RVM rendering to
% 'rvm_old'. Stephan
%
% 2018/01/19 Tor: Major update, adding canlab_load_ROI with region names
% and updating some atlases. Move towards all regions being precisely
% defined and usable as ROIs as well as display items.

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
  
    % -------------------------------------------------------------------
    % Subcortical regions using canlab_load_ROI
    % See canlab_load_ROI for more info and provenance of regions
    % -------------------------------------------------------------------
    
    case {'vmpfc' 'nacc' 'BST' ...
    'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
    'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
    'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
    'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'sc' 'ic' 'drn' 'mrn' ...
    'thalamus' 'thal' 'LGN' 'lgn' 'MGN' 'mgn' 'VPthal', 'VPLthal', 'VPL', 'intralaminar_thal', ...
    'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'ncf' 'vep'}
        
    [r, ~, default_color] = canlab_load_ROI(meth, 'noatlas');  % noatlas speeds things up! we don't need atlas 
    
    p = imageCluster('cluster',region2struct(r),'color',default_color,'alpha',.5);
    
    
    % -------------------------------------------------------------------
    % Surfaces
    % -------------------------------------------------------------------
    
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
    
    case 'hires surface left'
        pname = 'BigBrainSurfaceLeft.mat'; 
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);
     
        view(90,0);
        axis off;
        axis image;
        %         lightRestoreSingle(gca);
        material dull;
        
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


    case 'hires surface right'
        pname = 'BigBrainSurfaceRight.mat'; % from Glasser_et_al_2016_HCP
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);
     
        view(90,0);
        axis off;
        axis image;
        %         lightRestoreSingle(gca);
        material dull;
        
        
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
        
    case 'brainbottom'
        [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -20],'filename','scalped_single_subj_T1');
        set(p(1),'FaceColor',[.6 .4 .3]); colormap(gca, copper);material dull;axis off
        h = findobj('Type','Light'); delete(h); [az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
        set(p,'FaceAlpha',1)
        
    case {'cutaway', 'left_cutaway' 'right_cutaway' 'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5'}
        
        p = canlab_canonical_brain_surface_cutaways(meth, varargin{:});
        
        % Pre-2020: p = surface_cutaway('ycut_mm', -30);
        
    % -------------------------------------------------------------------
    % Other subcortical regions 
    % Some do not work well for display in canlab_load_ROI 
    % Some not in there because there are surface files, not regions
    % -------------------------------------------------------------------
            
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
        
    case 'amygdala'

        pname = 'surf_spm2_amy.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[0 0 .5]);
 case 'amygdala hires'

        pname = 'AMY.mat'; %anatomy toolbox

        p = add_surface(pname);
        set(p,'FaceColor',[0 0 .5]);
%     case 'thalamus'
%         %         P = which('carmack_thal_bstem.mat'); load(P)
%         %         p = imageCluster('cluster',thal,'color',[0 .8 .3],'alpha',.5);
%         pname = 'surf_spm2_thal.mat';
% 
%         p = add_surface(pname);
%         set(p,'FaceColor',[.9 .65 .5]);

    case {'hippocampus', 'hipp'}

        pname = 'surf_spm2_hipp.mat';

        p = add_surface(pname);
        set(p,'FaceColor',[.7 .4 .4]);
 
    case {'hippocampus hires'}

        pname = 'HC.mat'; %anatomy toolbox

        p = add_surface(pname);
        set(p,'FaceColor',[.7 .4 .4]);
%     case 'midbrain'
%         P = which('carmack_thal_bstem.mat');
%         load(P)
%         p = imageCluster('cluster',midbrain,'color',[.7 .3 0],'alpha',.5);
% 

    case {'cerebellum','cblm'}
        pname = 'surf_spm2_cblm.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[.8 .65 .8]);
     
    case {'CIT168'}
        pname ='CIT168.mat';
        p=add_surface(pname);
        set(p,'FaceColor',[.8 .65 .8]);
           
               
            
    % -------------------------------------------------------------------
    % Combinations of regions - batch
    % -------------------------------------------------------------------
    
    case {'BG', 'bg', 'basal ganglia'}
        
        names = {'caudate' 'put' 'GPe' 'GPi' 'VeP'};
        p = [];
        
        for i = 1:length(names)
            
            p = [p addbrain(names{i})];
            set(p(end), 'Tag', names{i});
            
        end

    case {'midbrain_group'}
        
        names = {'pag' 'sc' 'ic' 'drn' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn'};
        p = [];
        
        for i = 1:length(names)
            
            p = [p addbrain(names{i})];
            set(p(end), 'Tag', 'names');
            
        end
        
        view(140, 30)
        p = [p addbrain('brainstem')];
        set(p(end), 'FaceAlpha', .15);
        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        set(gca, 'ZLim', [-30 10])

    case 'limbic'
        p = [];
        myp = addbrain('amygdala');p = [p myp];
        myp = addbrain('hypothalamus');p = [p myp];
        myp = addbrain('hippocampus');p = [p myp];
        myp = addbrain('thalamus');p = [p myp];
        myp = addbrain('nacc');p = [p myp];
        myp = addbrain('hires surface left');p = [p myp];
         myp = addbrain('BG');p = [p myp];
         myp = addbrain('brainstem'); p = [p myp];
        set(p,'FaceAlpha',1);

        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
 case 'limbic hires'
        p = [];
        myp = addbrain('amygdala');p = [p myp];
        myp = addbrain('hippocampus');p = [p myp];
        myp = addbrain('hires left');p = [p myp];
        myp = addbrain('CIT168'); p = [p myp];
        myp = addbrain('brainstem'); p = [p myp];
        set(p,'FaceAlpha',1);

        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)

        
    case 'brainstem_group'
        
        names = {'mrn' 'pbn' 'lc' 'rvm' 'nts' 'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'pbn' 'ncf' }; % medulla only, midbrain separate
        % 'vep' 'PBP'
        p = addbrain('midbrain_group');
        
        for i = 1:length(names)
            
            p = [p addbrain(names{i})];
            set(p(end), 'Tag', 'names');
            
        end
        
        set(gca, 'ZLim', [-90 20])
        
    case {'thalamus_group'}
        
        names = {'lgn' 'mgn' 'VPthal', 'intralaminar_thal'};
        p = [];
        
        for i = 1:length(names)
            
            p = [p addbrain(names{i})];
            set(p(end), 'Tag', 'names');
            
        end
        
        view(140, 30)
        p = [p addbrain('brainstem')];
        p = [p addbrain('thalamus')];
        set(p(end), 'FaceAlpha', .15);
        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        set(gca, 'ZLim', [-30 10])
        
         
         
         % ---------------------------------
         
    case 'foursurfaces'
        p = run_foursurfaces;
        
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
    
    case {'limbic', 'BG', 'globus pallidus', 'gp', 'bg', 'basal ganglia', 'cutaway'}
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
surfh2 = copyobj(surfh, axh);
if iscolumn(surfh2), surfh2 = surfh2'; end
all_surf_handles = [all_surf_handles surfh2];

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
surfh2 = copyobj(surfh, axh);
if iscolumn(surfh2), surfh2 = surfh2'; end
all_surf_handles = [all_surf_handles surfh2];

view(90, 0);

lightRestoreSingle; axis image; axis off; lighting gouraud; material dull

end

