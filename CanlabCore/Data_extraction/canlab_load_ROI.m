function [r, obj, default_color, region_file, image_file] = canlab_load_ROI(region_name, varargin)
% Load a region by name (hand-picked from various atlases), for display or use as ROI in analysis
%
% - Easy to add regions from region objects or binary masks (.nii/.img)
% - Some regions best for display only; some good as ROIs as well
% - Inter-operates with addbrain.m to load regions for display
%
% :Usage:
% ::
%
%    handle = canlab_load_ROI(region_name,[optional arguments])
%
% Working options:
% -----------------------------------------------------------
% {'vmpfc' 'nacc' 'BST' ...
%     'cau' 'caudate' 'put' 'GP' 'GPe' 'GPi' 'VeP' ...
%     'thalamus' 'thal' 'cm' 'md' 'stn' 'habenula' 'mammillary' 'hypothalamus','hy','hythal' ...
%     'brainstem' 'midbrain' 'pag' 'PBP' 'sn' 'SNc' 'SNr' 'VTA' 'rn' ...
%     'pbn' 'lc' 'rvm' 'rvm_old' 'nts'}
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
% 'thalamus'
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
% 'pag'        Periaqueductal gray, hand-drawn (Tor Wager)
% 'PBP'        Parabrachial pigmented nuc.      % Pauli 2017 BioArxiv subcortical atlas
% 'sn'         Substantia Nigra; Keuken 2014
% 'SNc'        Substantia Nigra compacta        % Pauli 2017 BioArxiv subcortical atlas
% 'SNr'        Substantia Nigra reticularis     % Pauli 2017 BioArxiv subcortical atlas
% 'VTA'        Ventral tegmental area           % Pauli 2017 BioArxiv subcortical atlas
% 'rn'         Red nucleus; Keuken 2014
% 'pbn'        Parabrachial complex; very rough, hand-drawn for rough display (Tor)
% 'lc'         Locus coeruleus; Keren 2009, 2SD image
% 'rvm_old'    Hand-drawn rostral ventral medulla (Tor) in anatomical rvm
% 'rvm'        Rostral ventral medulla from Brooks et al. 2016(??)
% 'nts'        Nuc. tractus solitarius (rough; hand-drawn, Tor)
% 'olive'      Inferior olive; MISSING
% 'nrm'        Nuc. raphe magnus; MISSING
%
%
% Examples:
%
% [r, obj, default_color, region_file, image_file] = canlab_load_ROI('vmpfc');
% orthviews(r, 'color', {default_color});
%

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
    
else
    
    r = region();
    
end

% Load region object as r, convert to region object to standardize output
% -----------------------------------------------------------------------

if has_image
    
    % Load image file as obj, an fmri_data object
    
    obj = fmri_data(image_file);
    
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
        
    case 'thalamus'
        region_file = [];                                   % File with region object/clusters struct
        var_name = '';                                  % Variable name(s) of interest in file
        image_file = which('spm2_thal.img');   % Image file name with binary mask
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
        region_file = which('pag_cl.mat');              % File with region object/clusters struct
        var_name = 'pag';                               % Variable name(s) of interest in file
        image_file = which('spm5_pag.img');             % Image file name with binary mask
        default_color = [1 0 0];                        % default color for display
           
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
        region_file = which('pbn_cl.mat');
        var_name = 'pbn';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 .5 0];
        
    case {'lc', 'locus_coeruleus'}
        region_file = which('Keren_2009_LC_2SD_regions.mat');
        var_name = 'r';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 1 0];
        
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
        
    case 'nts'
        region_file = which('nts_cl.mat');
        var_name = 'nts';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [0 0 1];
        
    case {'olive', 'inferior olive'}
        region_file = [];
        var_name = '';                     % Variable name(s) of interest in file
        image_file = 'ROI_inf_olive.img';
        default_color = [.5 1 .5];
        
    case {'nrm', 'raphe magnus'}
        region_file = [];
        var_name = '';                     % Variable name(s) of interest in file
        image_file = 'ROI_raphe_magnus.img';
        default_color = [.3 0 1];
        
    otherwise
        error('Unknown region name.');
end

end % function
