function [r, atlas_obj, default_color, region_file, image_file] = canlab_load_ROI(region_name, varargin)
% Load a region by name (hand-picked from various atlases), for display or use as ROI in analysis
%
% - r is a region object, useful for display (e.g., addbrain.m and any region object method)
% - atlas_obj is an atlas object containing regions, mapped to a standard reference space (1 mm res)
%   this can be slow for some objects and is not needed for region display
%   (e.g., addblobs.m), so 'noatlas' will return an empty atlas_obj instead
%
% - Easy to add regions from region objects or binary masks (.nii/.img)
% - Some regions best for display only; some good as ROIs as well
% - Inter-operates with addbrain.m to load regions for display
%
%
% :Usage:
% ::
%
%    [r, atlas_obj, default_color, region_file, image_file] = canlab_load_ROI(region_name,[optional arguments])
%
% Working options:
% -----------------------------------------------------------
% {'vmpfc' 'nacc' 'BST' ...
%     'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
%     'thalamus' 'thal' 'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
%     'brainstem' 'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
%     'pbn' 'lc' 'rvm' 'rvm_old' 'nts' 'drn' 'mrn' 'sc' 'ic'}
%     
% Cortex -----------------------------------------------------------
% 'vmpfc'   Ventromedial prefrontal + posterior cing, midline; hand-drawn (Tor Wager)
%
% Forebrain (non-basal ganglia)
% -----------------------------------------------------------
% 'nacc'    Nucleus accumbens               % Pauli 2017 BioArxiv subcortical atlas
% 'hipp'    Hippocampus                     % MISSING/NEEDS UPDATE
% 'BST'     Bed nuc. of stria term/SLEA     % Pauli 2017 BioArxiv subcortical atlas
%
% Basal ganglia
% -----------------------------------------------------------
% 'caudate'  Caudate nucleus
% 'put'      Putamen            % MISSING
% 'GP'       Globus pallidus; Keuken 2014
% 'GPe'       Globus pallidus internal; Keuken 2014
% 'GPi'       Globus pallidus external; Keuken 2014
% 'VeP'        Ventral pallidum           % Pauli 2017 BioArxiv subcortical atlas
%
%  Thalamus, Diencephalon, Epithalamus
% -----------------------------------------------------------
% 'thalamus'                        Morel thalamus main body, Krauth 2010
% 'cm'      Centromedian thalamus   Morel thalamus atlas, Krauth 2010
% 'md'      Mediodorsal thalamus    Morel thalamus atlas, Krauth 2010
% 'stn'     subthalamic nucleus     Keuken 2014 (also options in Pauli, Morel)
% 'habenula'     Habenula           Pauli 2017 BioArxiv subcortical atlas (also in Morel)
% 'mammillary'   Mammillary bodies  Pauli 2017 BioArxiv subcortical atlas
% 'lgn'     Lateral geniculate nuc  Morel thalamus atlas, Krauth 2010
% 'mgn'     Medial geniculate nuc   Morel thalamus atlas, Krauth 2010
% 'VPthal'  Ventral posterior thal  Morel thalamus atlas, Krauth 2010
% 'intralaminar_thal' Intralaminar  Morel thalamus atlas, Krauth 2010
%
% Brainstem
% -----------------------------------------------------------
% 'brainstem'  Segmented and cleaned (Tor Wager) from SPM8 tissue probability maps
% 'midbrain'   Overall midbrain, from Carmack 2004
% 'pag'        Periaqueductal gray, hand-drawn (Tor Wager 2018, mask out aqueduct/Keuken2014)
% 'sc'         Superior colliculus, hand-drawn (Tor Wager 2018, mask out aqueduct/Keuken2014)
% 'ic'         Inferior colliculus, hand-drawn (Tor Wager 2018, mask out aqueduct/Keuken2014)
% 'drn'        Dorsal raphe nucleus, coords from Beliveau, 2015. Beliveau, Vincent, Claus Svarer, Vibe G. Frokjaer, Gitte M. Knudsen, Douglas N. Greve, and Patrick M. Fisher. 2015. ?Functional Connectivity of the Dorsal and Median Raphe Nuclei at Rest.? NeuroImage 116 (August): 187?95.
% 'mrn'        Median raphe nucleus, coords from Beliveau, 2015. Beliveau, Vincent, Claus Svarer, Vibe G. Frokjaer, Gitte M. Knudsen, Douglas N. Greve, and Patrick M. Fisher. 2015. ?Functional Connectivity of the Dorsal and Median Raphe Nuclei at Rest.? NeuroImage 116 (August): 187?95.
% 'PBP'        Parabrachial pigmented nuc.      % Pauli 2017 BioArxiv subcortical atlas
% 'sn'         Substantia Nigra; Keuken 2014   % Keuken, M. C., P-L Bazin, L. Crown, J. Hootsmans, A. Laufer, C. Müller-Axt, R. Sier, et al. 2014. ?Quantifying Inter-Individual Anatomical Variability in the Subcortex Using 7 T Structural MRI.? NeuroImage 94 (July): 40?46.
% 'SNc'        Substantia Nigra compacta        % Pauli 2017 BioArxiv subcortical atlas
% 'SNr'        Substantia Nigra reticularis     % Pauli 2017 BioArxiv subcortical atlas
% 'VTA'        Ventral tegmental area           % Pauli 2017 BioArxiv subcortical atlas
% 'rn'         Red nucleus; Keuken 2014
% 'pbn'        Parabrachial complex; Fairhurst, Merle, Katja Wiech, Paul Dunckley, and Irene Tracey. 2007. ?Anticipatory Brainstem Activity Predicts Neural Processing of Pain in Humans.? Pain 128 (1-2):101?10.
% 'lc'         Locus coeruleus; Keren 2009, 2SD image
% 'rvm_old'    Hand-drawn rostral ventral medulla (Tor) in anatomical rvm
% 'rvm'        Rostral ventral medulla from Brooks et al. 2016(??)
% 'nts'        Nuc. tractus solitarius (rough; hand-drawn, Tor)
% 'olive'      Inferior olive; MISSING
% 'nrm'        Nuc. raphe magnus; % Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% 'ncf'        Nuc. cuneiformis; Zambreanu, L., R. G. Wise, J. C. W. Brooks, G. D. Iannetti, and I. Tracey. 2005. ?A Role for the Brainstem in Central Sensitisation in Humans. Evidence from Functional Magnetic Resonance Imaging.? Pain 114 (3):397?407.
% 'ncs_B6_B8'  Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% 'nrp_B5'     Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% 'nuc_ambiguus' Sclocco, Roberta, Florian Beissner, Gaelle Desbordes, Jonathan R. Polimeni, Lawrence L. Wald, Norman W. Kettner, Jieun Kim, et al. 2016. ?Neuroimaging Brainstem Circuitry Supporting Cardiovagal Response to Pain: A Combined Heart Rate Variability/ultrahigh-Field (7 T) Functional Magnetic Resonance Imaging Study.? Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences 374 (2067). rsta.royalsocietypublishing.org. https://doi.org/10.1098/rsta.2015.0189.
% 'dmnx_nts'    Sclocco, Roberta, Florian Beissner, Gaelle Desbordes, Jonathan R. Polimeni, Lawrence L. Wald, Norman W. Kettner, Jieun Kim, et al. 2016. ?Neuroimaging Brainstem Circuitry Supporting Cardiovagal Response to Pain: A Combined Heart Rate Variability/ultrahigh-Field (7 T) Functional Magnetic Resonance Imaging Study.? Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences 374 (2067). rsta.royalsocietypublishing.org. https://doi.org/10.1098/rsta.2015.0189.
% 'vep'          % Pauli 2017 BioArxiv CIT168 subcortical atlas 
% 'medullary_raphe' Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.
% 'spinal_trigeminal' Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.
%
% Examples:
%
% [r, obj, default_color, region_file, image_file] = canlab_load_ROI('vmpfc');
% orthviews(r, 'color', {default_color});
%
% [r, myatlas] = canlab_load_ROI('habenula'); num_regions(myatlas); orthviews(myatlas)
% r = canlab_load_ROI('habenula', 'noatlas'); orthviews(r)

% Defaults and inputs
% -----------------------------------------------------------------------

doatlas = true;
reference_space_image = which('HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii');

if any(strcmp(varargin, 'noatlas')) || nargout < 2, doatlas = false; end

% Get names and info only first
[region_file, image_file, var_name, default_color] = get_region_info_from_registry(region_name);

has_region = ~isempty(region_file);
has_image = ~isempty(image_file);

% Load region object as r, convert to region object to standardize output
% -----------------------------------------------------------------------

if has_region
    
    if iscell(var_name)
        myregion = load(region_file, var_name{:});
        
        if isempty(myregion)
            disp('Cannot find region the variables you named in the file.');
            return
        end
        
        r = [];
        for i = 1:length(var_name)
            r = [r myregion.(var_name{i})];
        end
        
    else % not a cell
        
        myregion = load(region_file, var_name);
        
        if isempty(myregion)
            disp('Cannot find region the variables you named in the file.');
            return
        end
        
        r = myregion.(var_name);
        
    end

    if ~isa(r, 'region')
        r = cluster2region(r);
    end
    
    % Region names
    if strcmp(r(1).shorttitle, 'Region001')
        for i = 1:length(r)
            r(i).shorttitle = region_name;
            if iscell(var_name)
                r(i).shorttitle = var_name{1};
            else
                % r(i).shorttitle = var_name;
                r(i).shorttitle = region_name; % sometimes better...
                
            end
        end
    end
    
else
    
    r = region();
    
end

% Load region object as r, convert to region object to standardize output
% -----------------------------------------------------------------------

if has_image
    
    % Load image file as obj, an fmri_data object, to convert to region
    
    obj = fmri_data(image_file, 'noverbose');
    
else
    
    obj = fmri_data();
    
end


% If image but no region object yet, convert from image obj
% -----------------------------------------------------------------------

if has_image && ~has_region  % empty region
    
    r = region(obj);
    
end

% If region but no image, convert from region
% -----------------------------------------------------------------------

if ~has_image && has_region  % empty image
    
    try
        obj = region2fmri_data(r, obj);
        
    catch
        % cannot create obj if not in same space.  skip.
        
    end
end

% Re-map regions into standard reference space, 1 mm.  Needs image on path.
% -----------------------------------------------------------------------
atlas_obj = [];

if doatlas
    
    if isempty(reference_space_image)
        sprintf('To return atlas object, you need %s on your Matlab path. Skipping.\n',reference_space_image);
        return
    end
    
    atlas_obj = region2atlas(r, reference_space_image);
    
end

end % main function




function [region_file, image_file, var_name, default_color] = get_region_info_from_registry(region_name)

switch region_name
    
    % Cortex
    % -----------------------------------------------------------
    
    case {'vmpfc', 'vmPFC', 'VMPFC'}
        region_file = [];                               % File with region object/clusters struct
        var_name = '';                                  % Variable name(s) of interest in file
        image_file = which('VMPFC_display_mask.img');   % Image file name with binary mask
        default_color = [.7 .3 0];                      % default color for display
        
        
        % Forebrain
        % -----------------------------------------------------------
        
    case {'nucleus accumbens','nacc','nac'}
        %region_file = which('NucAccumb_clusters.mat');  % File with region object/clusters struct
        
        region_file = which('CIT168_atlas_regions.mat');  
        var_name = 'NAC';                                % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [0 .5 0];                       % default color for display
        
        
    case 'amygdala'
        region_file = which('amy_clusters.mat');         % File with region object/clusters struct
        var_name = 'amy';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [0 0 .5];                        % default color for display
        
    case {'hippocampus', 'hipp'} % ***************
        
        region_file = which('hipp_clusters.mat');  % File with region object/clusters struct
        var_name = '';                             % Variable name(s) of interest in file
        image_file = [];                           % Image file name with binary mask
        default_color = [.5 .6 .6];
                   
    case {'BST','BNST','SLEA'}

        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'BST_SLEA';                          % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [0 .5 1];                       % default color for display
      
      % Basal ganglia
        % -----------------------------------------------------------   
        
    case {'Cau', 'cau', 'caudate'}
        
        region_file = which('CIT168_atlas_regions.mat');  
        var_name = 'Cau';                                   % Variable name(s) of interest in file
        image_file = [];                                    % Image file name with binary mask
        default_color = [.5 .6 .6];
       
    case {'putamen', 'put'} 
        
        region_file = which('CIT168_atlas_regions.mat');  
        var_name = 'Put';                                   % Variable name(s) of interest in file
        image_file = [];                                    % Image file name with binary mask
        default_color = [.3 .7 .5];
        
    case {'globus pallidus', 'gp'}
        
        region_file = which('Keuken_2014_7T_regions.mat');  % File with region object/clusters struct
        var_name = {'GPe', 'GPi'};                          % Variable name(s) of interest in file
        image_file = [];                                    % Image file name with binary mask
        default_color = [.5 .6 .5];
        
    case {'GPe', 'gpe'}
        
        region_file = which('Keuken_2014_7T_regions.mat');  % File with region object/clusters struct
        var_name = 'GPe';                                   % Variable name(s) of interest in file
        image_file = [];                                    % Image file name with binary mask
        default_color = [.5 .6 .5];
        
    case {'GPi', 'gpi'}
        
        region_file = which('Keuken_2014_7T_regions.mat');  % File with region object/clusters struct
        var_name = 'GPi';                                   % Variable name(s) of interest in file
        image_file = [];                                    % Image file name with binary mask
        default_color = [.5 .75 .5];
        
        
        % Diencephalon and epithalamus
        % -----------------------------------------------------------
    case {'hypothalamus','hy','hythal'}

        %region_file = which('hy_clusters.mat');         % File with region object/clusters struct
        %var_name = 'hy';  
        
        region_file = which('CIT168_atlas_regions.mat');          
        var_name = 'Hythal';                                % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [1 1 0];                        % default color for display
        
    case {'HN', 'habenula'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'Haben';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .5 1];                       % default color for display
        
    case {'mamm', 'mammillary'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'Mamm_Nuc';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .5 1];                       % default color for display
        
        
        % Thalamus
        % -----------------------------------------------------------
        
    case {'thalamus', 'thal'}
        region_file = 'Morel_Thalamus_main_body_region_and_atlas_obj';  % File with region object/clusters struct
        var_name = 'thal_region';                                  % Variable name(s) of interest in file
        image_file = [];   % Image file name with binary mask.  old: which('spm2_thal.img');
        default_color = [.9 .65 .5];
        
        
    case {'md','mediodorsal'}
        
        region_file = which('Morel_thalamus_atlas_regions.mat');              % File with region object/clusters struct
        var_name = {'MD_group'};                                    % Variable name(s) of interest in file
        image_file = [];                                            % Image file name with binary mask
        default_color = [1 0 0];                                    % default color for display
        
    case {'cm','centromedian'}
        
        region_file = which('Morel_thalamus_atlas_regions.mat');              % File with region object/clusters struct
        var_name = {'L_CM', 'R_CM'};                               % Variable name(s) of interest in file
        image_file = [];             % Image file name with binary mask
        default_color = [1 0 0];                        % default color for display
        
    case {'stn', 'subthalamic nucleus'}
        region_file = which('Keuken_2014_7T_regions.mat');
        var_name = 'STN';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 0 .5];
        
    case {'LGN', 'lgn'}
        
        region_file = which('Morel_thalamus_atlas_regions.mat');    % File with region object/clusters struct
        var_name = {'LGN_group'};                                   % Variable name(s) of interest in file
        image_file = [];                                            % Image file name with binary mask
        default_color = [.5 0 1];                                   % default color for display
        
    case {'MGN', 'mgn'}
        
        region_file = which('Morel_thalamus_atlas_regions.mat');              % File with region object/clusters struct
        var_name = {'L_MGN' 'R_MGN'};                               % Variable name(s) of interest in file
        image_file = [];                                            % Image file name with binary mask
        default_color = [1 0 .5];                                   % default color for display
        
      case {'VPthal', 'VPLthal', 'VPL'}
        
        region_file = which('Morel_thalamus_atlas_regions.mat');    % File with region object/clusters struct
        var_name = {'VPL_group'};                                   % Variable name(s) of interest in file
        image_file = [];                                            % Image file name with binary mask
        default_color = [.3 .7 .3];                                   % default color for display      
        
       case {'intralaminar_thal'}
        
        region_file = which('Morel_thalamus_atlas_regions.mat');    % File with region object/clusters struct
        var_name = {'intralaminar_midline_group'};                                   % Variable name(s) of interest in file
        image_file = [];                                            % Image file name with binary mask
        default_color = [.8 .8 0];                                   % default color for display      
        
        
        
        
        % General brainstem
        % -----------------------------------------------------------
    case 'brainstem'
        
        region_file = [];                            % File with region object/clusters struct
        var_name = '';                               % Variable name(s) of interest in file
        image_file = which('spm8_brainstem.img');    % Image file name with binary mask
        default_color = [.5 .65 .4];                 % default color for display
        
        
        % Midbrain
        % -----------------------------------------------------------

    case 'midbrain'
        region_file = which('carmack_thal_bstem.mat');  % File with region object/clusters struct
        var_name = 'midbrain';                          % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.7 .3 0];                      % default color for display
        
    case 'pag'
        % 'ROI_pag.img' alternate
        region_file = which('pag_roi_2018_tor.mat');    % File with region object/clusters struct % old: 'pag_cl.mat'
        var_name = 'pag_regions';                       % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [1 0 0];                        % default color for display
        
    case 'sc'
        region_file = which('sup_inf_colliculus_roi_2018_tor.mat');  % File with region object/clusters struct
        var_name = 'sc_regions';                          % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.7 .5 .2];                      % default color for display
        
    case 'ic'
        region_file = which('sup_inf_colliculus_roi_2018_tor.mat');  % File with region object/clusters struct
        var_name = 'ic_regions';                          % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.4 .7 .2];                      % default color for display
        
    case 'ncf'
        region_file = which('coordinate_brainstem_rois_2018_tor.mat'); % old region_file = which('pbn_cl.mat');
        var_name = 'ncf_regions';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 .5 0];
        
    case 'drn'
        region_file = which('dorsal_median_raphe_roi_2018_tor.mat');  % File with region object/clusters struct
        var_name = 'drn_regions';                          % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .3 1];                      % default color for display
        
    case 'mrn'
        region_file = which('dorsal_median_raphe_roi_2018_tor.mat');  % File with region object/clusters struct
        var_name = 'mrn_regions';                          % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .3 1];                      % default color for display
        
    case {'PBP'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'PBP';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .5 1];                       % default color for display
        
    case {'sn', 'substantia nigra'}
        region_file = which('Keuken_2014_7T_regions.mat');
        var_name = 'SN';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [0 0 .5];
        
        
    case {'rn', 'red nucleus'}
        region_file = which('Keuken_2014_7T_regions.mat');
        var_name = 'RN';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 1 0];
        
        case {'SNc', 'snc'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'SNc';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .7 .3];                       % default color for display
        
        case {'SNr', 'snr'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'SNr';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .7 .5];                      % default color for display
        
        case {'VTA', 'vta'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'VTA';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .7 .2];                       % default color for display
        
        case {'vep', 'VeP', 'vpall'}
        
        region_file = which('CIT168_atlas_regions.mat');
        var_name = 'VeP';                               % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [.3 .5 1];                       % default color for display
        

        
        % Pons
        % -----------------------------------------------------------
        
    case 'pbn'
        region_file = which('coordinate_brainstem_rois_2018_tor.mat'); % old region_file = which('pbn_cl.mat');
        var_name = 'pbn_regions';                     % Variable name(s) of interest in file
         image_file = [];
        default_color = [1 .5 0];
        
    case {'lc', 'locus_coeruleus'}
        region_file = which('Keren_2009_LC_2SD_regions.mat');
        var_name = 'r';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 1 0];
        
    case {'nrp_B5'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'nrp_B5_regions';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [.3 .3 1];
        
        % Medulla
        % -----------------------------------------------------------
        
    case 'rvm_old' % legacy RVM mask
        region_file = which('rvm_cl.mat');  % File with region object/clusters struct
        var_name = 'r';                     % Variable name(s) of interest in file
        image_file = [];                        % Image file name with binary mask
        default_color = [1 .7 1];               % default color for display
        
    case {'rvm','rvm_brooks'} % RVM mask from Jon Brooks
        region_file = which('RVMmask_symm_2mm_Brooks.mat');
        var_name = 'r';                     % Variable name(s) of interest in file
        image_file = which('RVMmask_symm_2mm.nii');
        default_color = [1 .7 1];
        
%     case 'nts'
%         region_file = which('coordinate_brainstem_rois_2018_tor.mat');  % old which('nts_cl.mat');
%         var_name = 'dmnx_nts_regions';                     % Variable name(s) of interest in file
%         image_file = [];
%         default_color = [0 0 1];
        
    case {'olive', 'inferior olive'}
        region_file = [];
        var_name = '';                     % Variable name(s) of interest in file
        image_file = 'ROI_inf_olive.img';
        default_color = [.5 1 .5];
   
    case {'nrm', 'raphe magnus'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'nrm_regions';                     % Variable name(s) of interest in file
        image_file = [];                    % old: 'ROI_raphe_magnus.img';
        default_color = [.3 .3 1];
        
    case {'medullary_raphe'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'medullary_raphe_regions';                     % Variable name(s) of interest in file
        image_file = [];                    
        default_color = [.3 .3 1];
    
    case {'ncs_B6_B8'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'ncs_B6_B8_regions';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [.3 .3 1];
        
    case {'spinal_trigeminal'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'spinal_trigeminal_regions';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [.7 .3 .4];
        
    case {'nuc_ambiguus'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'nuc_ambiguus_regions';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [.8 .2 .2];
        
    case {'nts', 'dmnx_nts'}
        region_file = which('coordinate_brainstem_rois_2018_tor.mat');
        var_name = 'dmnx_nts_regions';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 .1 1];
        
        
    otherwise
        error('Unknown region name.');
end

end % function
