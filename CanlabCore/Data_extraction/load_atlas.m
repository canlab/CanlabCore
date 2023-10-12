function atlas_obj = load_atlas(atlas_file_name_or_keyword, varargin)
% Load one of a collection of atlases by keyword
%
% atlas_obj = load_atlas(varargin)
%
% List of keywords/atlases available:
% -------------------------------------------------------------------------
% 'canlab2018'                    'Combined atlas from other published atlases, whole brain. 
% 'canlab2018_2mm'                'Combined atlas resampled at 2 mm resolution'
% 'thalamus'                      'Thalamus_combined_atlas_object.mat'
% 'thalamus_detail', 'morel'      'Morel_thalamus_atlas_object.mat'
% 'cortex', 'glasser'             'Glasser2016HCP_atlas_object.mat'
% 'basal_ganglia', 'bg'           'Basal_ganglia_combined_atlas_object.mat'
% 'striatum', 'pauli_bg'          'Pauli2016_striatum_atlas_object.mat'
% 'brainstem'                     'brainstem_combined_atlas_object.mat'
% 'subcortical_rl', 'cit168'      'CIT168_MNI_subcortical_atlas_object.mat'
% 'brainnetome'                   'Brainnetome_atlas_object.mat'
% 'keuken'                        'Keuken_7T_atlas_object.mat'
% 'buckner'                       'buckner_networks_atlas_object.mat'
% 'cerebellum', 'suit'            'SUIT_Cerebellum_MNI_atlas_object.mat'
% 'shen'                          'Shen_atlas_object.mat'
% 'schaefer400'                   *Not saved as object yet* 'Schaefer2018Cortex_atlas_regions.mat' 
% 'yeo17networks'                 'Schaefer2018Cortex_17networks_atlas_object.mat'
% 'insula'                        'Faillenot_insular_atlas.mat'
% 'painpathways'                  'pain_pathways_atlas_obj.mat'
% 'painpathways_finegrained'      'pain_pathways_atlas_obj.mat'
% 'tian_3t_s1/2/3/4_fmriprep20/fsl6'      'Subcortical atlas at four different resolutions and two different reference spaces (e.g., tian_3t_s4_fmriprep20) )'
% 'delavega'                      'delaVega2017_neurosynth_atlas_object'
% 'cytoarchitecture_jubrain'      'Amunts2020_jubrain_cytoarchitecture.mat'
% 'jubrain/jubrain_fmriprep20/jubrain_fsl6' 'Julich Brain atlas in fmriprep 20.0.3 LTS (default) or fsl spaces
%
% More information and references to original publications are saved in
% each atlas object. This function is a shell to collect them in a central registry.
% New atlases can be created by passing a file name (e.g., .nii file) or an fmri_data object
% and labels into the atlas( ) constructor method.
%
% Examples:
% -------------------------------------------------------------------------
% atlas_obj = load_atlas('thalamus');
% atlas_obj = load_atlas('Thalamus_atlas_combined_Morel.mat');
%
%

docustom = 0;
verbose = 1;

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
    
    case {'thalamus'}
        savefile = which('Thalamus_combined_atlas_object.mat');
        varname = 'thalamus_atlas';
        
    case {'thalamus_detail', 'morel'}
        savefile = which('Morel_thalamus_atlas_object.mat');
        varname = 'atlas_obj';
        
    case {'cortex', 'glasser'}
        savefile = which('Glasser2016HCP_atlas_object.mat');
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
        
    case {'subcortical_rl', 'CIT168', 'cit168'}
        savefile = which('CIT168_MNI_subcortical_atlas_object.mat');
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
        
    case {'cerebellum', 'suit'}
        savefile = which('SUIT_Cerebellum_MNI_atlas_object.mat');
        varname = 'atlas_obj';
        
    case 'shen'
        savefile = which('Shen_atlas_object.mat');
        varname = 'atlas_obj';
        
        %          case 'schaefer400'
        %              savefile = which('Schaefer2018Cortex_atlas_regions.mat');
        %              varname = 'atlas_obj';
        
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
        
    case {'tian2020_subcortical_scale1', 'tian_3t_s1_fmriprep20_2mm'}
        savefile ='tian_3t_s1_fmriprep20_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian2020_subcortical_scale2', 'tian_3t_s2_fmriprep20_2mm'}
        savefile ='tian_3t_s2_fmriprep20_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian2020_subcortical_scale3', 'tian_3t_s3_fmriprep20_2mm'}
        savefile ='tian_3t_s3_fmriprep20_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian2020_subcortical_scale4', 'tian_3t_s4_fmriprep20_2mm'}
        savefile ='tian_3t_s4_fmriprep20_2mm_atlas_object.mat';
        varname = 'atlas_obj';
       
    case {'tian_3t_s1_fmriprep20'}
        savefile ='tian_3t_s1_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s2_fmriprep20'}
        savefile ='tian_3t_s2_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s3_fmriprep20'}
        savefile ='tian_3t_s3_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s4_fmriprep20'}
        savefile ='tian_3t_s4_fmriprep20_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s1_fsl6_2mm'}
        savefile ='tian_3t_s1_fsl6_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s2_fsl6_2mm'}
        savefile ='tian_3t_s2_fsl6_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s3_fsl6_2mm'}
        savefile ='tian_3t_s3_fsl6_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s4_fsl6_2mm'}
        savefile ='tian_3t_s4_fsl6_2mm_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s1_fsl6'}
        savefile ='tian_3t_s1_fsl6_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s2_fsl6'}
        savefile ='tian_3t_s2_fsl6_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s3_fsl6'}
        savefile ='tian_3t_s3_fsl6_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'tian_3t_s4_fsl6'}
        savefile ='tian_3t_s4_fsl6_atlas_object.mat';
        varname = 'atlas_obj';
        
    case {'jubrain','jubrain_fmriprep20'}
        savefile = 'julich_fmriprep20_atlas_object.mat';
        varname = 'juAtlas';

    case 'jubrain_fsl6'
        savefile = 'julich_fsl6_atlas_object.mat';
        varname = 'juAtlas';

    case 'delavega'
        savefile ='delaVega2017_neurosynth_atlas_object.mat';
        varname = 'atlas_obj';

    case {'cyto', 'jubrain', 'amunts2020'}
        warning('This atlas is deprecated. Please see ''jubrain'' instead.')
        savefile = which('Amunts2020_jubrain_cytoarchitecture.mat');
        varname = 'atlas_obj';
        
    otherwise % assume it's a file name
        savefile = which(atlas_file_name_or_keyword);
        varname = [];
        
end % switch

atlas_obj = load_atlas_from_file(savefile, varname, verbose);

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



