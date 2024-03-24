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
% % Regions (see below, and canlab_load_ROI.m)
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
% 'bigbrain' {'hires surface left', 'bigbrain left'}
% ['fsavg_left' or 'inflated left'], 'fsavg_right' or 'inflated right', 
% uses freesurfer inflated brain with Thomas Yeo group's RF_ANTs mapping from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
% 
% % Macro subcortical surfaces
% % -----------------------------------------------------------------------
% 'pauli_subcortical'    High-res combined surface from the Pauli CIT168 "reinforcement learning" atlas
% 'CIT168'               Lower-res combined surface from the Pauli CIT168 "reinforcement learning" atlas
% 'cerebellum','cblm'    surf_spm2_cblm,  an old but nice surface for rendering
% 'brainstem'            surf_spm2_brainstem, an old but nice surface for rendering
% 'suit brainstem'       Diedrichsen SUIT brainstem and cerebellum; 'suit_surface_brainstem_cerebellum.mat
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
%        a high-resolution cortical surface (from Caret segmentation)
%        of the Colin27 single-subject brain
%
%   **'hires left':**
%        hi-resolution left medial with cerebellum (Caret seg)
%
%   **'hires right':**
%        same, right hem
%
%   **'surface left':**
%        hi-resolution left cortical surface with cerebellum, defaults in medial view
%        
%        Pial surface with MSM alignment based on the Human Connectome Project
%        Reference: Glasser et al. 2016 Nature 
%        Surface template: L.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat'; 
%
%   **'surface right':**
%        same, right hem
%
%   **'hcp inflated left': same as the above but using inflated surfaces
%
%   **'hcp inflated right': same as the above but using inflated surfaces
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
%   **'left_cutaway'**      These surfaces use the Keuken 2014 7T MNI surface
%   **'right_cutaway' **
%   **'right_cutaway_x8' -- like right_cutaway but x=8
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
%   **'inflated surfaces'**
%        Left and right inflated cortical hemispheres, uses fsavg
%        Freesurfer with Yeo transformation to MNI
%
%   **'flat surfaces'**
%        Left and right inflated cortical hemispheres, uses fsavg
%        Freesurfer with Yeo transformation to MNI
%
%   **'insula surfaces'**
%        Left and right insula slabs and inflated surfaces
%        Zooms in on insula
%
%   **'multi_surface'**
%        Left and right cutaways (with basal ganglia/brainstem)
%        Left and right cortical surfaces (average brain)
%        Insula surfaces
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

% --------------------------------------------------------
% Create surfaces
% --------------------------------------------------------


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
            
            len = size(get(myp(i), 'FaceVertexCData'), 1);
            
            % If we have intensity data saved, use that
            % This is saved for isocaps
            usedprev = false;
            prevdata = [];
            
            try prevdata = get(myp(i), 'UserData'); catch, end % may not exist in all Matlab versions?
            
            if ~isempty(prevdata) && size(prevdata, 1) == len
                
                set(myp(i), 'FaceVertexCData', prevdata); % [.5 .5 .5]
                usedprev = true;
                
            else
                set(myp(i), 'FaceVertexCData', repmat(128, len, 1)); % [.5 .5 .5]
            end
            
        end
        p = myp;
        
        return
  
    % -------------------------------------------------------------------
    % Subcortical regions using canlab_load_ROI
    % See canlab_load_ROI for more info and provenance of regions
    % -------------------------------------------------------------------
    
    case {'BNST', 'bnst'}
        error('For Bed Nucleus of Stria Terminalis use BST')
        
    case {'vmpfc' 'nacc' 'BST' ...
    'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
    'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
    'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
    'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'sc' 'ic' 'drn' 'mrn' ...
    'thalamus' 'thal' 'LGN' 'lgn' 'MGN' 'mgn' 'VPthal', 'VPLthal', 'VPL', 'intralaminar_thal', ...
    'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'ncf' 'vep'}
        
    [r, ~, default_color] = canlab_load_ROI(meth, 'noatlas');  % noatlas speeds things up! we don't need atlas 
    
    p = imageCluster('cluster',region2struct(r),'color',default_color,'alpha',.5);
    
    if all(ishandle(p)), set(p, 'Tag', meth); end
        
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
    
    case {'hires surface left', 'bigbrain left'}
        
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


    case {'hires surface right', 'bigbrain right'}
        
        pname = 'BigBrainSurfaceRight.mat'; % from Amunts et al. Julich BigBrain
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5]);
     
        view(90,0);
        axis off;
        axis image;
        %         lightRestoreSingle(gca);
        material dull;
        
        
    case {'inflated right' 'fsavg_right'} % uses freesurfer inflated brain
            % with Thomas Yeo group's RF_ANTs mapping
            % from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
            pname = 'surf_freesurf_inflated_Right.mat';
             p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
     
        view(90, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'inflated left' 'fsavg_left'} 
        % uses freesurfer inflated brain with Thomas Yeo group's RF_ANTs mapping from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
        pname = 'surf_freesurf_inflated_Left.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
        
    case {'hcp inflated left'}
        % uses inflated midthickness image distributed with HCP_utils
        % package. By contrast 'inflated left' seems to use the 'very
        % inflated' image in HCP_utils, but I didn't create that so I'm not
        % changing the name.
        % - Bogdan        
        pname = 'S12000.L.inflated_MSMAll.32k_fsl_LR.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
        
    case {'hcp inflated right'}
        % uses inflated midthickness image distributed with HCP_utils
        % package. By contrast 'inflated left' seems to use the 'very
        % inflated' image in HCP_utils, but I didn't create that so I'm not
        % changing the name.
        % - Bogdan
        pname = 'S12000.R.inflated_MSMAll.32k_fsl_LR.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
    
    case {'hcp inflated'}
        % bilateral versions of the above
        pname = 'S12000.inflated_MSMAll.32k_fsl_LR.mat';
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'hcp sphere left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('S1200.L.sphere.32k_fs_LR.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'hcp sphere right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan     
        pname = which('S1200.R.sphere.32k_fs_LR.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'freesurfer sphere left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('fsavg_sphere_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'freesurfer sphere right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan     
        pname = which('fsavg_sphere_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'freesurfer inflated left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('fsavg_inflated_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'freesurfer inflated right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan     
        pname = which('fsavg_inflated_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'freesurfer white left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('fsavg_white_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'freesurfer white right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan     
        pname = which('fsavg_white_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
        
    case {'MNI152NLin2009cAsym white left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin2009cAsym_white_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin2009cAsym white right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin2009cAsym_white_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'MNI152NLin2009cAsym midthickness left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin2009cAsym_midthickness_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin2009cAsym midthickness right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin2009cAsym_midthickness_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'MNI152NLin2009cAsym pial left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin2009cAsym_pial_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin2009cAsym pial right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin2009cAsym_pial_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'MNI152NLin6Asym white left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_white_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin6Asym white right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_white_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'MNI152NLin6Asym midthickness left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_midthickness_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin6Asym midthickness right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_midthickness_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'MNI152NLin6Asym pial left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_pial_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin6Asym pial right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_pial_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case {'MNI152NLin6Asym sphere left'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_sphere_lh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;
            
    case {'MNI152NLin6Asym sphere right'}
        % Mean for use as a development tool, not intended for actual
        % dispaly
        % - Bogdan        
        pname = which('MNI152NLin6Asym_sphere_rh.mat');
        p = add_surface(pname);
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
        
        view(270, 0);
        axis off;
        axis image;
        lightRestoreSingle(gca);
        material dull;

    case 'bigbrain'
        
        pname = 'BigBrainSurfaceRight.mat'; 
        p = add_surface(pname);
        
        pname = 'BigBrainSurfaceLeft.mat'; 
        p = [p add_surface(pname)];
        
        set(p,'FaceColor',[.5 .5 .5], 'FaceAlpha', 1);
     
        view(135, 30);
        axis off;
        axis image;
        material dull;
        lightRestoreSingle;
        
        
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
        
    case {'cutaway', 'left_cutaway' 'right_cutaway' 'right_cutaway_x8' 'left_insula_slab' 'right_insula_slab' 'accumbens_slab' 'coronal_slabs' 'coronal_slabs_4' 'coronal_slabs_5'}
        
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
        % Lower-res surface from the Pauli CIT168 "reinforcement learning" atlas
        pname ='CIT168.mat';
        p=add_surface(pname);
        set(p,'FaceColor',[.8 .65 .8]);
           
    case {'pauli_subcortical'}
        % High-res surface from the Pauli CIT168 "reinforcement learning" atlas
        pname ='pauli_subcortical.mat';
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
        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca);
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

        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca);
        
 case 'limbic hires'
        p = [];
        myp = addbrain('amygdala');p = [p myp];
        myp = addbrain('hippocampus');p = [p myp];
        myp = addbrain('hires left');p = [p myp];
        myp = addbrain('CIT168'); p = [p myp];
        myp = addbrain('brainstem'); p = [p myp];
        set(p,'FaceAlpha',1);

        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca);

        
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
        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca);
        set(gca, 'ZLim', [-30 10])
        
    case 'inflated surfaces'
        clf;
        subplot(1, 2, 1);
        p = addbrain('fsavg_left');
        subplot(1, 2, 2);
        p = [p addbrain('fsavg_right')];
        
    case 'flat surfaces'
        clf;
        subplot(1, 2, 1);
        p = addbrain('flat left');
        subplot(1, 2, 2);
        p = [p addbrain('flat right')];
      
    case 'insula surfaces'
        clf;
        subplot(2, 2, 1);
        p = addbrain('left_insula_slab');
        axis off
        set(gca, 'YLim', [-40 40], 'ZLim', [-25 25]);
        
        subplot(2, 2, 2);
        p = [p addbrain('right_insula_slab')];
        axis off
        set(gca, 'YLim', [-40 40], 'ZLim', [-25 25]);
        
        subplot(2, 2, 3);
        p = [p addbrain('fsavg_left')];
        camzoom(1.5)
        set(gca, 'YLim', [-40 40], 'ZLim', [-25 25]);
        
        ax = subplot(2, 2, 4);
        p = [p addbrain('fsavg_right')];
        axes(ax)
        camzoom(1.5)
        set(gca, 'YLim', [-40 40], 'ZLim', [-25 25]);  % not working for some reason - don't know why!! camzoom works
        disp('set(gca, ''YLim'', [-40 40], ''ZLim'', [-25 25]);');
        
        
        % four surfaces, inflated and flat
        
    case 'multi_surface'
        
        create_figure('multi_surface', 2, 2);
        subplot(2, 2, 1);
        p = addbrain('right_cutaway');
        camzoom(1.3)
        
        subplot(2, 2, 2);
        p = [p addbrain('left_cutaway')];
        camzoom(1.3)
        
        subplot(2, 2, 3);
        p = [p addbrain('surface left')];
        view(270,0);
        lightRestoreSingle
        camzoom(1.3)
        
        subplot(2, 2, 4);
        p = [p addbrain('surface right')];
        lightRestoreSingle
        camzoom(1.3)

        set(p, 'FaceAlpha', 1);
        drawnow
        
        create_figure('multi_surface2', 2, 2);
        p = [p addbrain('insula surfaces')];
        
        set(p, 'FaceAlpha', 1);

        
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

% --------------------------------------------------------
% Set Tag
% --------------------------------------------------------

switch meth
    
    case {'limbic', 'BG', 'globus pallidus', 'gp', 'bg', 'basal ganglia', 'cutaway', 'brainstem_group', 'limbic hires', 'left_cutaway', 'right_cutaway'}
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

