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

% Cortex -----------------------------------------------------------
% 'vmpfc'   Ventromedial prefrontal + posterior cing, midline; hand-drawn (Tor Wager)
%
% Forebrain (non-basal ganglia)
% -----------------------------------------------------------
% 'nacc'    Nucleus accumbens   % Drawn by Tor around Knutson coords I think (provenance uncertain). Best for region display only.
% hipp'     Hippocampus         % MISSING/NEEDS UPDATE
%
% Basal ganglia
% -----------------------------------------------------------
% 'caudate'  Caudate nucleus
% 'put'      Putamen            % MISSING
% 'GP'       Globus pallidus; Keuken 2014
% 'GPe'       Globus pallidus internal; Keuken 2014
% 'GPi'       Globus pallidus external; Keuken 2014
%
%  Thalamus and Diencephalon
% -----------------------------------------------------------
% 'thalamus'
% 'cm'      Centromedian thalamus   Hand-drawn by Tor; Best for region display only.
% 'md'      Mediodorsal thalamus    Hand-drawn by Tor; Best for region display only.
% 'stn'     subthalamic nucleus     Keuken 2014
%
% Brainstem
% -----------------------------------------------------------
% 'brainstem'  Segmented and cleaned (Tor Wager) from SPM8 tissue probability maps
% 'midbrain'   Overall midbrain, from Carmack 2004
% 'pag'        Periaqueductal gray, hand-drawn (Tor Wager)
% 'sn'         Substantia Nigra; Keuken 2014
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
    else
        myregion = load(region_file, var_name);
    end
    
    r = myregion.(var_name);
    
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
        region_file = which('NucAccumb_clusters.mat');  % File with region object/clusters struct
        var_name = 'cl';                                % Variable name(s) of interest in file
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
        
      % Basal ganglia
        % -----------------------------------------------------------   
        
    case 'caudate' % ***************
        %P = which('Tal_Cau.img'); %which('ICBM_caudate.img');   %('Tal_Cau.img'); %
        %[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
        %if findstr(P,'Tal_Cau.img'), str(49:58) = '[.5 .6 .6]'; delete(p(1)); p = p(2);,eval(str),end
        %         P = which('carmack_more_clusters.mat'); load(P)
        %         p = imageCluster('cluster',cau,'color',[.5 .6 .6],'alpha',.5);
        
        %         pname = 'surf_spm2_caudate.mat';
        %
        %         p = add_surface(pname);
        %         set(p,'FaceColor',[.5 .6 .6]);
        
                
    case {'putamen', 'put'} % ***********************
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
        %         pname = 'spm_surf_putamen_luke.mat';
        %          p = add_surface(pname);
        %         set(p,'FaceColor',[.3 .7 .5]);
        
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
        
        
        % Diencephalon
        % -----------------------------------------------------------
    case {'hypothalamus','hy','hythal'}
        
        region_file = which('hy_clusters.mat');         % File with region object/clusters struct
        var_name = 'hy';                                % Variable name(s) of interest in file
        image_file = [];                                % Image file name with binary mask
        default_color = [1 1 0];                        % default color for display
        
        % Thalamus
        % -----------------------------------------------------------
        
    case 'thalamus'
        region_file = [];                                   % File with region object/clusters struct
        var_name = '';                                  % Variable name(s) of interest in file
        image_file = which('spm2_thal.img');   % Image file name with binary mask
        default_color = [.9 .65 .5];
        
        
    case {'md','mediodorsal'}
        
        region_file = which('thal_brainstem_approx_working.mat');              % File with region object/clusters struct
        var_name = 'MD';                               % Variable name(s) of interest in file
        image_file = [];             % Image file name with binary mask
        default_color = [1 0 0];                        % default color for display
        
    case {'cm','centromedian'}
        
        region_file = which('thal_brainstem_approx_working.mat');              % File with region object/clusters struct
        var_name = 'CM';                               % Variable name(s) of interest in file
        image_file = [];             % Image file name with binary mask
        default_color = [1 0 0];                        % default color for display
        
    case {'stn', 'subthalamic nucleus'}
        region_file = which('Keuken_2014_7T_regions.mat');
        var_name = 'STN';                     % Variable name(s) of interest in file
        image_file = [];
        default_color = [1 0 .5];
        
        
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
        
        % Pons
        % -----------------------------------------------------------
        
    case 'pbn'
        region_file = which('pbn_cl.mat');
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
        var_name = 'r';                     % Variable name(s) of interest in file
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
