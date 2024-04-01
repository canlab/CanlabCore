function atlas_obj = load_atlas(atlas_file_name_or_keyword, varargin)
% Load one of a collection of atlases by keyword
%
% atlas_obj = load_atlas(varargin)
%
% List of keywords/atlases available:
% -------------------------------------------------------------------------
% 'canlab2024[_fine|_coarse][_fmriprep20|_fsl6][_1mm|_2mm]
%                                 'Combined atlas from other published, whole brain. Available in a fine or coarse (default) parcellation in
%                                 MNI152NLin2009cAsym (aka fmriprep) space (default) or MNI152NLin6Asym (aka fsl) space in 1mm or 2mm (default)
%                                 resolutions. Additional parcellations available with downsample_parcellation(). Refer to github README for 
%                                 details. Assembled dynamically when called due to distribution restrictions of Bianciardi subatlas. Development 
%                                 is ongoing (Date: 03/07/2024)'
% 'opencanlab2024[_fine|_coarse][_fmriprep20|_fsl6][_1mm|2mm_]
%                                 'Variation on canlab2024 that only uses regions with open usage and distribution licenses. No client side 
%                                 assembly needed, but has fewer brainstem areas, and some have synthetic probabilities. Used as the starting
%                                 point for dynamic assembly of canlab2024 though and the two have much in common. Defauts = {coarse,fmriprep20,2mm}
% 'canlab2023[_fine|_coarse][_fmriprep20|_fsl6][_1mm|_2mm]
%                                 'Combined atlas from other published, whole brain. Available in a fine or coarse (default) parcellation in
%                                 MNI152NLin2009cAsym (aka fmriprep) space (default) or MNI152NLin6Asym (aka fsl) space in 1mm or 2mm (default)
%                                 resolutions. Additional parcellations available with downsample_parcellation(). Refer to github README for 
%                                 details. Development is frozen, so it may be more stable than CANLab2024 (Date: 3/07/2024)'
% 'canlab2018[_2mm]'              'Combined atlas from other published atlases, whole brain. (Deprecated in favor of canlab2023)' 
% 'desikan_killiany[_fsl6|_fmriprep20]'
%                                 'Desikan-Killiany cortical gyral/sulcal labeling from Freesurfer (2006). Gyri include pial surface and lateral banks. Default: fmriprep20'
% 'dkt[_fsl6|_fmriprep20]         'Kline and Tourville's update to the desikan killiany atlas projected to fsl or fmriprep spaces. Default: fmriprep20'
% 'destrieux[_fsl6|_fmriprep20]'  'Destrieux cortical gyral/sulcal labeling from Freesurfer (2009). Gyral/sulcal division is based on curvature values. Default: fmriprep20'
% 'thalamus'                      'Thalamus_combined_atlas_object.mat'
% 'thalamus_detail', 'morel[_fsl6|_fmriprep20]',
%                                 'Morel_thalamus_atlas_object.mat in MNI152NLin6Asym (fsl) space (default) or MNI152NLin2009cAsym (fmriprep) space. 
%                                 (Both in MasksPrivate)'
% 'iglesias_thal[_fmriprep20|_fsl6]'
%                                 'Iglesias/Freesurfer thalamic nuclear parcellation in fmriprep20 (default) or fsl6 space. A bit more coarse than morel, 
%                                  but open license, more accurate boundaries and probablistic'
% 'iglesias_hypothal[_fmriprep20|_fsl6]
%                                 'Billot/Iglesias/Freesurfer hypothalamic segmentation in fmriprep20 (default) or fsl6 space.
% 'cortex', 'glasser'
%                                 'Glasser 2016 multimodal cortical parcellation volumetric projection using nearest neighbor interpolation from 
%                                  surface (deprecated)'
% 'glasser_[fmriprep20|fsl6]'     'Glasser 2016 multimodal cortical parcellation volumetric projection using registration fusion to two surface 
%                                  templates using 2 studies (N=241/89)
% 'basal_ganglia', 'bg'           'Basal_ganglia_combined_atlas_object.mat'
% 'striatum', 'pauli_bg'          'Pauli2016_striatum_atlas_object.mat'
% 'brainstem'                     'brainstem_combined_atlas_object.mat'
% 'subcortical_rl','cit168'       'CIT168 MNI152Nlin2009cAsym subcortical atlas v1.0.0 (deprecated)'
% 'cit168_[fmriprep20|fsl6]'      'CIT168 v1.1.0 subcortical atlas in fmriprep20 or fsl6 space'
% 'cit168_amygdala_[fmriprep20|fsl6]
%                                 'CIT168 v1.0.3 amygdalar nuclear parcellation in fmriprep20 or fsl6 space'
% 'brainnetome'                   'Brainnetome_atlas_object.mat'
% 'keuken'                        'Keuken_7T_atlas_object.mat'
% 'buckner'                       'buckner_networks_atlas_object.mat'
% 'cerebellum[_fsl6|_fmriprep20]', 'suit[_fsl6|_fmriprep20]'
%                                 'Diedrichsen cerrebellar atlas in MNI152NLin6Asym space (aka fsl, default) or MNI152NLin2009cAsym space (aka fmriprep).'
% 'shen[_fmriprep20|_fsl6]'       'Shen_atlas_object.mat, in MNIColin27v1998 space (default), MNI152NLin6Asym (fsl) and MNI152NLin2009cAsym (fmriprep 20.2.3) spaces'
% 'schaefer400'                   *Not saved as object yet* 'Schaefer2018Cortex_atlas_regions.mat' 
% 'yeo17networks'                 'Schaefer2018Cortex_17networks_atlas_object.mat'
% 'insula'                        'Faillenot_insular_atlas.mat'
% 'painpathways'                  'pain_pathways_atlas_obj.mat'
% 'painpathways_finegrained'      'pain_pathways_atlas_obj.mat'
% 'tian_3t_[fmriprep20|fsl6]'      
%                                 'Subcortical atlas at four different resolutions and two different reference spaces. Use atlas/get_coarser_parcellation to select low resolution versions.'
% 'delavega'                      'delaVega2017_neurosynth_atlas_object'
% 'julich_[fmriprep20|fsl6]'      'Histological Julich Brain atlas in fmriprep 20.2.3 LTS (default) or fsl spaces'
% 'bianciardi[_fmriprep20|_fsl6][_2mm]   
%                                 'Bianciardi brainstem atlas in fmriprep 20.2.3 LTS space (default) or fsl spaces at 1mm (default) or 2mm sampling resolution'
% 'cartmell_NAc[_fmriprep20|_fsl6] 'NAc Core/Shell probablistic atlas'
% 'harvard_aan[_fmriprep20|_fsl6] 'Harvard ascending arousal network atlas version 2.0. A generalization of the reticular activating system to various other brainstem nuclei besides the 
%                                   midbrain reticular formation. Based on histology, immunihistochemistry, and DWI tractography.'
% 'limbic_brainstem_atlas[_fmriprep20|fsl_6]
%                                 'Levinson Bari Limbic Brainstem Atlas. Includes VTA, dorsal raphe, locus coereleus, nucleus tractus solitaris and PAG. Probablistic with an open usage license.'
%
% More information and references to original publications are saved in
% each atlas object. This function is a shell to collect them in a central registry.
% New atlases can be created by passing a file name (e.g., .nii file) or an fmri_data object
% and labels into the atlas( ) constructor method.
%
% Examples:
% -------------------------------------------------------------------------
% atlas_obj = load_atlas('thalamus');
% atlas_obj = load_atlas('Thalamus_atlas_combined_Morel.mat');morel
%
%

docustom = 0;
verbose = 1;
docreate = 0;

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'noverbose', verbose = 0;
                
            case 'verbose', verbose = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

switch lower(atlas_file_name_or_keyword)
    
    case {'thalamus', 'morel'}
        warning('This atlas is deprecated. Please specify thalamus_[fsl6|fmriprep20] instead.');
        savefile = which('Thalamus_combined_atlas_object.mat');
        varname = 'thalamus_atlas';
        
    case {'thalamus_detail', 'morel_fsl6'}
        savefile = which('Morel_thalamus_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'morel_fmriprep20'}
        savefile = which('Morel_MNI152NLin2009cAsym_atlas_object.mat');
        varname = 'atlas_obj';

    case {'iglesias_thal_fmriprep20', 'iglesias_thal'}
        savefile = which('iglesias_HCP278_ST76_PG264_MNI152NLin2009cAsym_atlas_object.mat');
        varname = 'atlas_obj';

    case {'iglesias_thal_fsl6'}
        savefile = which('iglesias_HCP278_ST76_PG264_MNI152NLin6Asym_atlas_object.mat');
        varname = 'atlas_obj';

    case {'iglesias_hypothal_fmriprep20', 'iglesias_hypothal'}
        savefile = which('iglesias_hypothal_hcp278_MNI152NLin2009cAsym_atlas_object.mat');
        varname = 'atlas_obj';

    case {'iglesias_hypothal_fsl6'}
        savefile = which('iglesias_hypothal_hcp278_MNI152NLin6Asym_atlas_object.mat');
        varname = 'atlas_obj';

    case {'cortex', 'glasser'}
        warning('This version of the Glasser atlas is deprecated. Please invoke glasser_[fmriprep20|fsl6]');
        savefile = which('Glasser2016HCP_atlas_object.mat');
        varname = 'atlas_obj';

    case 'glasser_fmriprep20'
        savefile = which('glasser_MNI152NLin2009cAsym_atlas_object.mat');
        varname = 'atlas_obj';

    case 'glasser_fsl6'
        savefile = which('glasser_MNI152NLin6Asym_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'basal_ganglia', 'bg'}
        savefile = which('Basal_ganglia_combined_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'brainstem'}
        savefile = which('brainstem_combined_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'striatum', 'pauli_bg'}
        savefile = which('Pauli2016_striatum_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'subcortical_rl', 'cit168'}
        warning('This version (v1.0.0) is deprecated. Please use v1.1.0 by invoking cit168_fmriprep20 or cit168_fsl6.')
        savefile = which('CIT168_MNI_subcortical_atlas_object.mat');
        varname = 'atlas_obj';

    case {'cit168_fmriprep20'}
        savefile = which('CIT168_MNI152NLin2009cAsym_subcortical_v1.1.0_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'cit168_fsl6'}
        savefile = which('CIT168_MNI152NLin6Asym_subcortical_v1.1.0_atlas_object.mat');
        varname = 'atlas_obj';

    case {'cit168_amygdala_fmriprep20', 'cit168_amygdala'}
        savefile = which('CIT168_MNI152NLin2009cAsym_amygdala_v1.0.3_atlas_object.mat');
        varname = 'atlas_obj';

    case {'cit168_amygdala_fsl6'}
        savefile = which('CIT168_MNI152NLin6Asym_amygdala_v1.0.3_atlas_object.mat');
        varname = 'atlas_obj';

    case {'brainnetome'}
        savefile = which('Brainnetome_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'keuken'}
        savefile = which('Keuken_7T_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'buckner'}
        %savefile = 'Buckner1000FC_2011_cortex_atlas_object.mat';
        savefile = 'buckner_networks_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'cerebellum', 'suit', 'cerebellum_fsl6', 'suit_fsl6'}
        savefile = which('SUIT_Cerebellum_MNI_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'cerebellum_fmriprep20', 'suit_fmriprep20'}
        savefile = which('SUIT_Cerebellum_MNI152NLin2009cAsym_atlas_object.mat');
        varname = 'atlas_obj';

    case 'shen'
        warning('This is Shen in MNIColin27v1998 space, a subject specific space of the original paper. Consider shen_[fmriprep20|fsl6] instead.')
        savefile = which('Shen_atlas_object.mat');
        varname = 'atlas_obj';
        
        %          case 'schaefer400'
        %              savefile = which('Schaefer2018Cortex_atlas_regions.mat');
        %              varname = 'atlas_obj';
    case 'shen_fmriprep20'
        savefile = which('Shen_MNI152NLin2009cAsym_atlas_object.mat');
        varname = 'atlas_obj';
        
    case 'shen_fsl6'
        savefile = which('Shen_MNI152NLin6Asym_atlas_object.mat');
        varname = 'atlas_obj';
        
    case 'yeo17networks'
        savefile = which('Schaefer2018Cortex_17networks_atlas_object.mat');
        varname = 'atlas_obj';
        
    case 'canlab2018'
        savefile = 'CANlab_combined_atlas_object_2018.mat';
        varname = 'atlas_obj';
        
    case 'canlab2018_2mm'
        savefile = 'CANlab_combined_atlas_object_2018_2mm.mat';
        varname = 'atlas_obj';

    case 'insula'
        savefile = 'Faillenot_insular_atlas.mat';
        varname = 'atlas_obj';
        
    case 'painpathways'
        savefile = 'pain_pathways_atlas_obj.mat';
        varname = 'pain_pathways';
        
    case 'painpathways_finegrained'
        savefile = 'pain_pathways_atlas_obj.mat';
        varname = 'pain_pathways_finegrained';
        
    case 'kragel2019pag'
        savefile ='Kragel2019PAG_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t', 'tian_3t_fmriprep20'}
        savefile ='tian_3t_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_fsl6'}
        savefile ='tian_3t_fsl6_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'julich','julich_fmriprep20'}
        savefile = 'julich_fmriprep20_atlas_object.mat';
        varname = 'juAtlas';

    case 'julich_fsl6'
        savefile = 'julich_fsl6_atlas_object.mat';
        varname = 'juAtlas';

    case 'delavega'
        savefile ='delaVega2017_neurosynth_atlas_object.mat';
        varname = 'atlas_obj';

    case {'cartmell_nac', 'cartmell_nac_fmriprep20'}
        savefile='NAcCoreShell_MNI152NLin2009cAsym_atlas_object.mat';
        varname='atlas_obj';
        
    case 'cartmell_nac_fsl6'
        savefile='NAcCoreShell_MNI152NLin6Asym_atlas_object.mat';
        varname='atlas_obj';
        
    case {'bianciardi', 'bianciardi_fmriprep20'}
        savefile='bianciardi_MNI152NLin2009cAsym_atlas_object.mat';
        varname = 'bianciaAtlas';
        docreate = true;
        create_atlas = @(x1)bianciardi_create_atlas_obj('MNI152NLin2009cAsym');
        
    case {'bianciardi_fsl6'}
        savefile='bianciardi_MNI152NLin6Asym_atlas_object.mat';
        varname = 'bianciaAtlas';
        docreate = true;
        create_atlas = @(x1)bianciardi_create_atlas_obj('MNI152NLin6Asym');

    case {'bianciardi_2mm', 'bianciardi_fmriprep20_2mm'}
        savefile='bianciardi_MNI152NLin2009cAsym_2mm_atlas_object.mat';
        varname = 'bianciaAtlas';
        docreate = true;
        create_atlas = @(x1)bianciardi_create_atlas_obj('MNI152NLin2009cAsym_2mm');

    case {'bianciardi_fsl6_2mm'}
        savefile='bianciardi_MNI152NLin6Asym_2mm_atlas_object.mat';
        varname = 'bianciaAtlas';
        docreate = true;
        create_atlas = @(x1)bianciardi_create_atlas_obj('MNI152NLin6Asym_2mm');
        
    case {'canlab2023_coarse_fmriprep20_2mm', 'canlab2023_coarse_fmriprep20','canlab2023_coarse_2mm', ...
            'canlab2023_fmriprep20_2mm', 'canlab2023_coarse', 'canlab2023_fmriprep20', 'canlab2023_2mm', ...
            'canlab2023'}
        savefile='CANLab2023_MNI152NLin2009cAsym_coarse_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin2009cAsym','coarse',2);

    case {'canlab2023_coarse_fsl6_2mm', 'canlab2023_coarse_fsl6', 'canlab2023_fsl6_2mm', 'canlab2023_fsl6'}
        savefile='CANLab2023_MNI152NLin6Asym_coarse_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin6Asym','coarse',2);

    case {'canlab2023_coarse_fmriprep20_1mm', 'canlab2023_coarse_1mm', 'canlab2023_fmriprep20_1mm', 'canlab2023_1mm'}
        savefile='CANLab2023_MNI152NLin2009cAsym_coarse_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin2009cAsym','coarse',1);

    case {'canlab2023_coarse_fsl6_1mm', 'canlab2023_fsl6_1mm'}
        savefile='CANLab2023_MNI152NLin6Asym_coarse_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin6Asym','coarse',1);

    case {'canlab2023_fine_fmriprep20_2mm', 'canlab2023_fine_fmriprep20','canlab2023_fine_2mm', 'canlab2023_fine'}
        savefile='CANLab2023_MNI152NLin2009cAsym_fine_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin2009cAsym','fine',2);

    case {'canlab2023_fine_fsl6_2mm', 'canlab2023_fine_fsl6'}
        savefile='CANLab2023_MNI152NLin6Asym_fine_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin6Asym','fine',2);
        
    case {'canlab2023_fine_fmriprep20_1mm', 'canlab2023_fine_1mm'}
        savefile='CANLab2023_MNI152NLin2009cAsym_fine_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin2009cAsym','fine',1);

    case {'canlab2023_fine_fsl6_1mm'}
        savefile='CANLab2023_MNI152NLin6Asym_fine_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2023_atlas('MNI152NLin6Asym','fine',1);

    case {'opencanlab2024_fine_fmriprep20_1mm', 'opencanlab2024_fine_1mm'}
        savefile='openCANLab2024_MNI152NLin2009cAsym.mat';
        varname = 'canlab';

    case {'opencanlab2024_fine_fmriprep20_2mm', 'opencanlab2024_fine_fmriprep20', 'opencanlab2024_fine'}
        savefile='openCANLab2024_MNI152NLin2009cAsym_2mm.mat';
        varname = 'canlab_2mm';

    case {'opencanlab2024_fine_fsl6_1mm'}
        savefile='openCANLab2024_MNI152NLin6Asym.mat';
        varname = 'canlab';

    case {'opencanlab2024_fine_fsl6_2mm', 'opencanlab2024_fine_fsl6'}
        savefile='openCANLab2024_MNI152NLin6Asym_2mm.mat';
        varname = 'canlab_2mm';

    case {'opencanlab2024_coarse_fmriprep20_1mm', 'opencanlab2024_coarse_1mm', 'opencanlab2024_1mm'}
        savefile='openCANLab2024_MNI152NLin2009cAsym_coarse.mat';
        varname = 'canlab_coarse';

    case {'opencanlab2024_coarse_fmriprep20_2mm', 'opencanlab2024_coarse_fmriprep20', 'opencanlab2024_coarse_2mm', 'opencanlab2024_fmriprep20_2mm', ...
            'opencanlab2024_coarse', 'opencanlab2024_fmriprep20', 'opencanlab2024_2mm', 'opencanlab2024'}
        savefile='openCANLab2024_MNI152NLin2009cAsym_coarse_2mm.mat';
        varname = 'canlab_coarse_2mm';

    case {'opencanlab2024_coarse_fsl6_1mm', 'opencanlab2024_fsl6_1mm'}
        savefile='openCANLab2024_MNI152NLin6Asym_coarse.mat';
        varname = 'canlab_coarse';

    case {'opencanlab2024_coarse_fsl6_2mm', 'opencanlab2024_fsl6_2mm', 'opencanlab2024_coarse_fsl6', 'opencanlab2024_fsl6'}
        savefile='openCANLab2024_MNI152NLin6Asym_coarse_2mm.mat';
        varname = 'canlab_coarse_2mm';

    case {'canlab2024_coarse_fmriprep20_2mm', 'canlab2024_coarse_fmriprep20','canlab2024_coarse_2mm', ...
            'canlab2024_fmriprep20_2mm', 'canlab2024_coarse', 'canlab2024_fmriprep20', 'canlab2024_2mm', ...
            'canlab2024'}
        savefile='CANLab2024_MNI152NLin2009cAsym_coarse_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin2009cAsym','coarse',2);

    case {'canlab2024_coarse_fsl6_2mm', 'canlab2024_coarse_fsl6', 'canlab2024_fsl6_2mm', 'canlab2024_fsl6'}
        savefile='CANLab2024_MNI152NLin6Asym_coarse_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin6Asym','coarse',2);

    case {'canlab2024_coarse_fmriprep20_1mm', 'canlab2024_coarse_1mm', 'canlab2024_fmriprep20_1mm', 'canlab2024_1mm'}
        savefile='CANLab2024_MNI152NLin2009cAsym_coarse_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin2009cAsym','coarse',1);

    case {'canlab2024_coarse_fsl6_1mm', 'canlab2024_fsl6_1mm'}
        savefile='CANLab2024_MNI152NLin6Asym_coarse_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin6Asym','coarse',1);

    case {'canlab2024_fine_fmriprep20_2mm', 'canlab2024_fine_fmriprep20','canlab2024_fine_2mm', 'canlab2024_fine'}
        savefile='CANLab2024_MNI152NLin2009cAsym_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin2009cAsym','fine',2);

    case {'canlab2024_fine_fsl6_2mm', 'canlab2024_fine_fsl6'}
        savefile='CANLab2024_MNI152NLin6Asym_fine_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin6Asym','fine',2);
        
    case {'canlab2024_fine_fmriprep20_1mm', 'canlab2024_fine_1mm'}
        savefile='CANLab2024_MNI152NLin2009cAsym_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin2009cAsym','fine',1);

    case {'canlab2024_fine_fsl6_1mm'}
        savefile='CANLab2024_MNI152NLin6Asym_fine_atlas_object.mat';
        varname = 'atlas_obj';
        docreate = true;
        create_atlas = @(x1)create_CANLab2024_atlas('MNI152NLin6Asym','fine',1);

    case {'harvard_aan', 'harvard_aan_fmriprep20'}
        savefile='harvard_aan_v2_MNI152NLin2009cAsym_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'harvard_aan_fsl6'}
        savefile='harvard_aan_v2_MNI152NLin6Asym_atlas_object.mat';
        varname = 'atlas_obj';

    case {'limbic_brainstem_atlas_fmriprep20', 'limbic_brainstem_atlas'}
        savefile='levinson_bari_limbic_brainstem_atlas_MNI152NLin2009cAsym_atlas_object.mat';
        varname = 'atlas_obj';

    case {'limbic_brainstem_atlas_fsl6'}
        savefile='levinson_bari_limbic_brainstem_atlas_MNI152NLin6Asym_atlas_object.mat';
        varname = 'atlas_obj';

    case {'desikan_killiany', 'desikan_killiany_fmriprep20'}
        savefile='desikan_killiany_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';

    case {'desikan_killiany_fsl6'}
        savefile='desikan_killiany_fsl6_atlas_object.mat';
        varname = 'atlas_obj';

    case {'dkt', 'dkt_fmriprep20'}
        savefile='dkt_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';

    case {'dkt_fsl6'}
        savefile='dkt_fsl6_atlas_object.mat';
        varname = 'atlas_obj';

    case {'destrieux', 'destrieux_fmriprep20'}
        savefile='destrieux_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';

    case {'destrieux_fsl6'}
        savefile='destrieux_fsl6_atlas_object.mat';
        varname = 'atlas_obj';


    otherwise % assume it's a file name
        savefile = which(atlas_file_name_or_keyword);
        varname = [];
        
end % switch

if ~docreate

    atlas_obj = load_atlas_from_file(savefile, varname, verbose);

else
    % has creation script, so let's run it if either is true: (1) atlas
    % file is missing, (2) atlas is out of date
    if isempty(dir(which(savefile))) 
        % atlas is missing
        create_atlas(0);
        atlas_obj = load_atlas_from_file(savefile, varname, verbose);
    else
        latest = which(strrep(savefile,'.mat','.latest'));
        if exist(latest, 'file') ~= 2
            create_atlas(0);
            latest = which(strrep(savefile,'.mat','.latest'));
        end
        fid = fopen(latest);
        latest = char(fread(fid,inf)');
        fclose(fid);

        atlas_obj = load_atlas_from_file(savefile, varname, verbose);
        fnames = fieldnames(atlas_obj.additional_info);
        if ~any(ismember('hash',fnames)) || ~strcmp(latest,atlas_obj.additional_info.hash)
            % evaluation is rounded to the nearest second to deal with
            % floating point issues. If there's an atlas 1s newer availble
            % then it means we rebuild the atlas to get it.
            if verbose
                fprintf('Updating atlas...\n');
            end
            create_atlas(0);
            
            atlas_obj = load_atlas_from_file(savefile, varname, verbose);
        end
    end
end


end % function


function atlas_obj = load_atlas_from_file(savefile, varname, verbose)

atlas_obj = [];

if isempty(savefile)
    fprintf('Cannot find atlas file: %s\n', savefile);
    return
end

if verbose
    fprintf('Loading atlas: %s\n', savefile);
end

if isempty(varname)
    
    tmp = load(savefile);
    vn = fieldnames(tmp);
    atlas_obj = tmp.(vn{1});
    
else
    
    atlas_obj = load(savefile, varname);
    
    if ~isfield(atlas_obj, varname), fprintf('Cannot find variable: %s\n', varname); return, end
    
    atlas_obj = atlas_obj.(varname);
    
end


end % subfunction



