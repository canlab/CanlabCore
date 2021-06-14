function [p,str] = cluster_surf(varargin)
% Surface plot of clusters on a standard brain
%
%    WARNING: SURFACE_CUTAWAY, TOR_3D, AND CLUSTER_SURF ARE DEPRECATED.
%    THEY USE AN OLDER STYLE OF SURFACE RENDERING THAT IS VERY SLOW. THEY
%    STILL WORK, BUT OBJECT-ORIENTED CANLAB TOOLS HAVE SHIFTED TO USING
%    ADDBRAIN.M COMBINED WITH RENDER_ON_SURFACE() METHOD, WHICH USES
%    MATLAB'S ISOCOLORS FOR DRAMATICALLY FASTER COLOR RENDERING.
%
% :Usage:
% ::
%
%    [suface_handle,colorchangestring] = cluster_surf(varargin)
%
% :Inputs:
%
%   **CLUSTERS:**
%        clusters structures, as created in tor_extract_rois.m
%
%   **COLORS:**
%        cell array of colors for each cluster: {[1 0 0] [0 1 0] [0 0 1]}
%       if number of colors specified is greater than number of clusters
%       structures entered, n+1 and n+2 colors are overlap of 2 and overlap
%       of all clusters, respectively.
%
%   **SURFACE MAT FILE:**
%        file name of mat file containing brain surface vertices and faces
%       as created with isosurface.
%
%   **SPECIAL SURFACE KEYWORDS:**
%        special string: 'bg' 'hipp' (hcmp,thal,amy)
%        number of mm to plot from surface (mmdeep)
%
%          - Special keywords for sets of surfaces are:
%            left, right, bg, limbic, cerebellum, brainstem
%
%          - Other keywords: 'left' 'right' 'amygdala' 'thalamus' 'hippocampus' '
%            'midbrain'  'caudate'   'globus pallidus'  'putamen'  'nucleus accumbens'
%             'hypothalamus' 'cerebellum'
%
%   **EXISTING SURFACE HANDLE(S):**
%        handles for surface patches, created,
%        e.g., with addbrain.m.  This lets you be very flexible in the
%        surfaces you image onto.
%
%   **'colorscale':**
%        This scales colors by Z-scores of voxels if used
%          - Uses input color, unless also used with 'heatmap'
%          - Z scores should be in ROW vector
%          - use with 'normalize' to scale Z-scores between -1 and 1
%          - will also create transparent effects, mixing
%            blob color with existing surface color in linear
%            proportion to Z-scores
%
%   **'heatmap':**
%        Map Z-scores to surface colors
%          - Used WITH or instead of 'colorscale'
%          - Blobs can have a range of colors
%          - Use with REFERENCE RANGE option below to control scale
%          - solid colors entered as input will be ignored
%          - use with 'colormaps' option below to be flexible
%            in which color maps you apply.
%          - if 'colorscale' is also used, will produce transparent blobs.
%
%   **REFERENCE RANGE:**
%        reference Z-scores range, [zmin_act zmax_act
%        zmax_negact zmin_negact], e.g., [0 5 -5 0], use only
%        with 'heatmap' option
%        to get refZ from clusters, try:
%        ::
%
%            clZ = cat(2,clusters.Z);
%            refZ = [min(clZ(clZ > 0)) max(clZ) min(clZ(clZ < 0)) min(clZ)];
%
%   **'colormaps':**
%        - followed by custom [colors x 3] matrices for positive colors
%          and negative colors.
%        - matlab can create some: e.g., colormap summer, jet, etc.
%          others can be created with colormap_tor.m
%
%        color [0 1 1] (cyan) is reserved for the overlap color btwn cluster sets.
%
% :Examples:
% ::
% % ------------------------------------------------------------
% % ------------------------------------------------------------
%
%    P = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_single_subj_T1_gray.mat';
%    cluster_surf(tcl,acl,P,10,{[0 1 0] [1 0 0]},'colorscale','heatmap')
%
%    or P = h (surface handle) to use current surface in figure, and refZ
%    cluster_surf(tcl,acl,h,[3 5 -5 -3],10,{[0 1 0] [1 0 0]},'colorscale','heatmap')
%
% % A complete example visualizing a network on several surfaces
% % ------------------------------------------------------------
% % Load atlas
% atl = load_atlas('yeo17networks');
% 
% % Select a network
% r = atlas2region(select_atlas_subset(atl, 1));
% 
% % Check what it looks like on slices
% figure; montage(r);
% 
% % Check what it looks like using surface() method
% % and render_on_surface()
% 
% figure; surface(r);
% 
% %% Run cluster_surf on several standard left-hemisphere surfaces
% % ------------------------------------------------------------
% create_figure('Left lateral surface', 2, 2);
% cluster_surf(r, 2, 'colors', {[1 1 0]}, 'hires left');
% view(275, 10); lightRestoreSingle;
% title('MNI Colin single-subj, addbrain(...''hires left'')');
% 
% subplot(2, 2, 2);
% % Inflated surfaces: These use the pre-calculated mappers from MNI to an
% % inflated freesurfer cortical surface for each hemisphere:
% % which('lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat')
% % which('rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat')
% try
%     cluster_surf(r, 2, 'colors', {[1 1 0]}, 'fsavg_left');
%     view(275, 10); lightRestoreSingle;
%     title('Freesurfer average with Wu 2017 MNI->surf mapping');
% catch
%     disp('cannot find required surface file');
% end
% 
% subplot(2, 2, 3);
% % Pial surface from Glasser_et_al_2016 Nature, HCP data
% % 'L.pial_MSMAll_2_d41_WRN_DeDrift.32k.mat';
% try
%     surf_han = addbrain('surface left');
%     cluster_surf(r, 2, 'colors', {[1 1 0]}, surf_han);
%     view(275, 10); lightRestoreSingle;
%     set(surf_han, 'FaceAlpha', 1);
%     title('HCP pial surface, addbrain(...''surface left'')');
% catch
%     title('cannot find required surface file');
% end
% 
% subplot(2, 2, 4);
% % BigBrain surface, Amunts et al.
% % 'BigBrainSurfaceLeft.mat' in addbrain
% try
%     surf_han = addbrain('hires surface left');
%     cluster_surf(r, 2, 'colors', {[1 1 0]}, surf_han);
%     view(275, 10); lightRestoreSingle;
%     set(surf_han, 'FaceAlpha', 1);
%     title('MNI Colin single-subj, addbrain(...''hires surface left'')');
% catch
%     title('cannot find required surface file');
% end
% % ------------------------------------------------------------
% % ------------------------------------------------------------
% :More examples:
% ::
% % ------------------------------------------------------------
%
%    cluster_surf(cl,2,'heatmap');     % brain surface.  vertices colored @2 mm
%    cluster_surf(cl,2,'bg','heatmap');    % heatmap on basal ganglia
%    cluster_surf(cl,5,'left','heatmap');  % heatmap on left surface @5 mm
%    cluster_surf(cl,2,'right','heatmap');
%
%    % A multi-color, multi-threshold display on the cerebellum
%    colors = {[1 1 0] [1 .5 0] [1 .3 .3]};
%    tor_fig;
%    sh = cluster_surf(cl{3},colors(3),5,'cerebellum');
%    cluster_surf(cl{2},colors(2),5,sh);
%    cluster_surf(cl{1},colors(1),5,sh);
%
%    % Custom colormaps:
%    create_figure('Brain Surface'); cluster_surf(cl, 2, 'heatmap','left');
%
%    poscm = colormap_tor([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
%    negcm = colormap_tor([0 0 1], [0 .3 1]);  % light blue to dark blue
%    create_figure('Brain Surface'); cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, 'left');
%
%    % Single-color transparent map (green):
%    cluster_surf(cl, 2, {[0 1 0]}, 'colorscale', p3(2), 'normalize');
%
% :See Also: addbrain, img2surf.m, surface() methods for objects,
% cluster_cutaways, render_on_surface()
%
%    WARNING: SURFACE_CUTAWAY, TOR_3D, AND CLUSTER_SURF 
%    USE AN OLDER STYLE OF SURFACE RENDERING THAT IS VERY SLOW. THEY
%    STILL WORK, BUT OBJECT-ORIENTED CANLAB TOOLS HAVE SHIFTED TO USING
%    ADDBRAIN.M COMBINED WITH RENDER_ON_SURFACE() METHOD, WHICH USES
%    MATLAB'S ISOCOLORS FOR DRAMATICALLY FASTER COLOR RENDERING.
%    CLUSTER_SURF IS STILL USEFUL FOR INFLATED SURFACES, AND FOR RENDERING
%    MULTIPLE SETS OF REGIONS ON THE SAME BRAIN IN DIFFERENT COLORS, AND
%    WITH FLEXIBLE CONTROL OVER THE DISTANCE TO SURFACE.

% ..
%    Programmers' Notes
%    Created by Tor, a long time ago
%    updated Sept 2015 to keep up with matlab graphics and handle some weird
%    stuff with processing inputs.
%      - Figures are now scalars and we need to check for those first.
%      - Also changed default surface and colormap
%
%    Jan 2020: Updated to add 'noverbose' option, minor cosmetic code
%    cleanup
% ..

disp('Cluster_surf uses a type of surface rendering with cluster_surf that is deprecated')
disp('It still functions, but the new method based on isocolors is dramatically faster')

% -----------------------------------------------------------------------
%    set up input arguments and defaults
% -----------------------------------------------------------------------
mmdeep = 10;
cscale = 0;
heatm = 0;
viewdeg = [135 30];
cl = [];

P = which('surf_spm2_brain.mat');
P = which('surf_spm2_brain_1mm.mat');

% default color maps
% These match fmridisplay:
poscm = colormap_tor([1 .5 0], [1 1 0]);
negcm = colormap_tor([0 0 1], [0 1 1]);

actcolors = [];  % used with heatmap
donormalize = 0; % used with colorscale
adjust_var = [];

doverbose = true;
if any(strcmp(varargin,'noverbose')), doverbose = false; end % do first

% -----------------------------------------------------------------------
%    optional inputs
% -----------------------------------------------------------------------
clind = 1;
for i = 1:length(varargin)
    
    if isempty(varargin{i})
        % ignore it
        
    elseif isstruct(varargin{i}) || isa(varargin{i}, 'region')
        cl{clind} = varargin{i}; clind = clind+1;
        
    elseif iscell(varargin{i}), mycolors = varargin{i};
        
    elseif isstr(varargin{i})
        
        if strcmp(varargin{i},'noverbose'), doverbose = false;
            
        elseif strcmp(varargin{i},'colorscale'), cscale = 1;
            
        elseif strcmp(varargin{i},'normalize'), donormalize = 1;
            
        elseif strcmp(varargin{i},'heatmap'), heatm = 1;
            
        elseif strcmp(varargin{i},'colormaps')
            if doverbose, disp('Using custom color maps.'); end
            poscm = varargin{i + 1}; varargin{i + 1} = [];
            negcm = varargin{i + 2};  varargin{i + 2} = [];
            
        elseif  strcmp(varargin{i},'left')
            P = which('surf_spm2_left.mat'); %which('surf_single_subj_grayL.mat');
            viewdeg = [90 0];
            
        elseif  strcmp(varargin{i},'right')
            P = which('surf_spm2_right.mat'); %which('surf_single_subj_grayR.mat');
            viewdeg = [270 0];
            
        elseif  strcmp(varargin{i},'hires left')
            P = which('surf_spm2_brain_left.mat'); %which('surf_single_subj_grayL.mat');
            viewdeg = [90 0];
            
        elseif  strcmp(varargin{i},'hires right')
            P = which('surf_spm2_brain_right.mat'); %which('surf_single_subj_grayR.mat');
            viewdeg = [270 0];
        elseif strcmp(varargin{i},'fsavg_right') % uses freesurfer inflated brain
            % with Thomas Yeo group's RF_ANTs mapping
            % from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
            P = which('surf_freesurf_inflated_Right.mat');
            viewdeg = [270 0];
            adjust_var = 'fsavg_right'; % varargin for getVertexColors
        elseif strcmp(varargin{i},'fsavg_left') % uses freesurfer inflated brain
            P = which('surf_freesurf_inflated_Left.mat');
            viewdeg = [90 0];
            adjust_var = 'fsavg_left'; % varargin for getVertexColors
        else P = varargin{i};
        end
        
    elseif isnumeric(varargin{i}) && isscalar(varargin{i})
        % it's a number, mm deep to render (any face within mmdeep mm of a significant voxel will
        % be colored)
        mmdeep = varargin{i};
        
    elseif any(ishandle(varargin{i})) && ~all(ishandle(varargin{i}))
        % Note: some weird things are happening as fig handles are class
        % double, and not recognized as figs, but are figure handles.
        % hopefully this will work with strcmp()
        % do nothing, but fig handles may be passed in inadvertently?
        % sometimes passing in integers is interpreted as fig handles...
        % For example, passing in 0 as a double is interpreted as a root
        % graphics object. Not a good design feature of new graphics in
        % Matlab...
        
        % Do nothing here. Could be all numeric input.
        
    elseif any(ishandle(varargin{i})) && any(strcmp(get(varargin{i}, 'Type'), 'figure'))
        % all are handles, one is a figure
        
        disp('You passed in a figure handle.')
        
    elseif  all(ishandle(varargin{i})) && all(isa(varargin{i}, 'matlab.graphics.primitive.Patch'))  %% all(~isa(varargin{i}, 'matlab.ui.Figure'))
        % all are handles, and patches
        % OLD: Pre-2014:  all(ishandle(varargin{i})) && all(varargin{i} ~= round(varargin{i}))
        if doverbose, disp('Found surface patch handles - plotting on existing surfaces.'); end
        P = varargin{i}; % handle(s) for existing surface
        
        % get rid of later calls to other surfaces
        wh = find(strcmp(varargin,'left') | strcmp(varargin,'right') | strcmp(varargin,'hires left') | strcmp(varargin,'hires right'));
        for jj = 1:length(wh), varargin{wh(jj)} = []; end
        
    elseif any(ishandle(varargin{i})) && ~all(isa(varargin{i}, 'matlab.graphics.primitive.Patch'))
        % OLD:  Pre-2014: any(ishandle(varargin{i})) && all(varargin{i} ~= round(varargin{i}))
        % Note: some weird things are happening as fig handles are class
        % double, and not recognized as figs, but are figure handles.
        % hopefully this will work with strcmp()
        if doverbose
            disp('Found existing surface patch handles, but some are invalid.  Check handles.');
            error('Exiting')
        end
        
    elseif length(varargin{i}) > 1   % it's a vector
        refZ = varargin{i};
        
    else    % no idea
        warning('cluster_surf: Unknown input');
        
    end
    
end

if ~exist('cl', 'var') || isempty(cl)
    if doverbose, disp('cluster_surf.m: No clusters to plot.  Try addbrain for brain surfaces with no activation.'); end
    p = []; str = [];
    return
end

if ~exist('mycolors', 'var')
    mycolors = scn_standard_colors(length(cl));
end

if isempty(P)
    disp(['Cannot find: ' P]);
    P = spm_get(1,'*mat','Choose brain surface file');
end

if doverbose
    
    disp('cluster_surf')
    disp('___________________________________________')
    fprintf('\t%3.0f cluster structures entered\n',length(cl))
    disp('  Colors are:')
    
    for i = 1:length(mycolors)
        disp([' ' num2str(mycolors{i})])
    end
    
end

if length(mycolors) > length(cl)
    if doverbose, disp([' overlap color is ' num2str(mycolors{length(cl)+1})]), end
    ovlc = ['[' num2str(mycolors{length(cl)+1}) ']'];
else
    ovlc = '[0 1 1]';
end

if length(mycolors) > length(cl)+1 && length(cl) > 2
    if doverbose, disp([' all overlap color is ' num2str(mycolors{length(cl)+2})]), end
    aovlc = ['[' num2str(mycolors{length(cl)+2}) ']'];
else
    aovlc = '[1 1 1]';
end

%disp([' Surface stored in: ' P])
if doverbose
    fprintf(' Building XYZ coord list\n');
end

% -------------------------------------------------------------------------
% * build xyz list
%
% also get cscale values for each coordinate, and alphascale values if both
% heatmap and colorscale options are entered
% -------------------------------------------------------------------------
for i = 1:length(cl) % each cell is a whole vector of cl, not a single region
    
    xyz{i} = cat(2,cl{i}.XYZmm)';
    
    if cscale || heatm
        for j = 1:length(cl{i})
            if size(cl{i}(j).Z,1) > size(cl{i}(j).Z,2)
                cl{i}(j).Z = cl{i}(j).Z';
            end
        end
        Z{i} = cat(2,cl{i}.Z)';
        
        % order voxels from lowest to highest, so that peak colors
        % appear because they are plotted last
        tmp = [xyz{i} Z{i}];
        tmp = sortrows(tmp,4);
        xyz{i} = tmp(:,1:3);
        Z{i} = tmp(:,4);1 ./ mad(abs(Z{i}));
        
        if cscale
            if donormalize
                Z{i} = Z{i} ./ max(Z{i});
            end
            
            if heatm
                % treat colorscale as alpha scaling to add transparent blobs
                % (preserve existing surface; good for isosurface objects)
                
                Za = abs(Z{i});
                a = prctile(Za, 85);  % midpoint of Z{i} defines transparency 0.5
                b = 1./mad(Za); % multiplier for Z for sigmoid
                
                sZ = 1 ./ (1 + exp(-b*(Za-a)));
                
                % fix, if all constant
                sZ(isnan(sZ)) = 1;
                
                alphascale{i} = sZ;
                
                % this is further adusted based on radius
                alphascale{i} = 5 * alphascale{i} ./ mmdeep^3; % should be 3?
                %alphascale
                
                
                
            end
            
        else
            % if heat map only, set mycolor{1} = [1 1 1]
            mycolors{1} = [1 1 1];
        end
    end
end

% ------------------------------------------------------------
% for heatmap option: get actcolors
% -------------------------------------------------------------
if heatm
    if doverbose, fprintf(' Getting heat-mapped colors\n'); end
    if exist('refZ') == 1
        actcolors = get_actcolors(Z, refZ, poscm, negcm);
    else
        actcolors = get_actcolors(Z, [], poscm, negcm);
    end
end

% rearrange Z to map the negative peak later
if exist('Z', 'var')
    for i = 1:numel(Z)
        neg_idx = Z{i}<0;
        xyz{i}(neg_idx,:)=flipud(xyz{i}(neg_idx,:));
        Z{i}(neg_idx)=flipud(Z{i}(neg_idx));
        actcolors{i}(neg_idx,:)=flipud(actcolors{i}(neg_idx,:));
    end
end

% -------------------------------------------------------------------------
% * build function call
% -------------------------------------------------------------------------
if doverbose, fprintf(' Building color change function call\n'); end

if length(cl) > 2 && exist('aovlc') == 1
    str = ['[c,alld] = getVertexColors(xyz{1},p,mycolors{1},[.5 .5 .5],' num2str(mmdeep) ',''ovlcolor'',' ovlc ',''allcolor'',' aovlc];
else
    str = ['[c,alld] = getVertexColors(xyz{1},p,mycolors{1},[.5 .5 .5],' num2str(mmdeep) ',''ovlcolor'',' ovlc];
end

if heatm
    str = [str ',''colorscale'',actcolors{1}'];
    if cscale
        % treat colorscale as alpha scaling to add transparent blobs
        str = [str ',''alphascale'',alphascale{1}'];
    end
elseif cscale
    % cscale alone - treat cscale as color-mapping index for single
    % colors in actcolors
    str = [str ',''colorscale'',Z{1}'];
end


for i = 2:length(cl)
    str = [str ',''vert'',xyz{' num2str(i) '},mycolors{' num2str(i) '}'];
    if heatm
        str = [str ',''colorscale'',actcolors{' num2str(i) '}'];
        if cscale
            % treat colorscale as alpha scaling to add transparent blobs
            str = [str ',''alphascale'',alphascale{' num2str(i) '}'];
        end
    elseif cscale
        % cscale alone - treat cscale as color-mapping index for single
        % colors in actcolors
        str = [str ',''colorscale'',Z{' num2str(i) '}'];
    end
end

if ~isempty(adjust_var)
    str = [str ',''' adjust_var ''''];
end

if ~doverbose
    str = [str ', ''noverbose'''];
end

str = [str ');'];

%p = P(end);
%co = get(p, 'FaceVertexCData');
if exist('alphascale','var')
    alphascale{1} = alphascale{1} * 12;
end
%eval(str)
%set(p, 'FaceVertexCData', co);

% -------------------------------------------------------------------------
% * run brain surface
% -------------------------------------------------------------------------
if ishandle(P)      % no input file, use existing handle
    if doverbose
        fprintf(' Using existing surface image\n');
        fprintf(' Running color change.\n');
    end
    for i = 1:length(P)
        p = P(i);
        if doverbose, disp([' eval: ' str]), end
        eval(str)
    end
else
    % we have either an input file or a special string ('bg')
    if doverbose, fprintf(' Loading surface image\n'); end
    [dtmp,ftmp,etmp]=fileparts(P);
    
    if strcmp(etmp,'.mat')
        
        load(P);
        
        %%figure
        p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
            'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
        lighting gouraud;camlight right
        axis image;
        lightRestoreSingle(gca);
        %myLight = camlight(0,0);set(myLight,'Tag','myLight');
        %set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView
        
        view(viewdeg(1),viewdeg(2));
        drawnow
        
        
        % -------------------------------------------------------------------------
        % * run color change
        % -------------------------------------------------------------------------
        if doverbose
            fprintf(' Running color change.\n');
            disp([' eval: ' str])
        end
        eval(str);
        
        
        % this for subcortex stuff
    elseif strcmp(P,'bg')
        p = [];
        myp = addbrain('caudate');p = [p myp];
        run_colorchange(myp,str,xyz,mycolors);
        
        myp = addbrain('globus pallidus');p = [p myp];
        run_colorchange(myp,str,xyz,mycolors);
        
        myp = addbrain('putamen');p = [p myp];
        run_colorchange(myp,str,xyz,mycolors);
        
        set(myp,'FaceAlpha',1);
        
        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
        
    elseif strcmp(P,'limbic')
        p = [];
        
        myp = addbrain('amygdala');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('hypothalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('hippocampus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('thalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('nucleus accumbens');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        
        myp = addbrain('left');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        set(myp,'FaceAlpha',1);
        
        axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
    elseif strcmp(P,'brainstem')
        p = [];
        
        myp = addbrain('thalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('hypothalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('brainstem');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        %myp = addbrain('caudate');p = [p myp];
        %run_colorchange(myp,str,xyz, mycolors, actcolors);
        
        
        view(90,10); axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
    elseif strcmp(P,'subcortex')
        p = [];
        
        myp = addbrain('thalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('hypothalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('brainstem');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('amygdala');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('hippocampus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('thalamus');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('nucleus accumbens');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('caudate');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        myp = addbrain('putamen');p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        %         myp = addbrain('globus pallidus');p = [p myp];  %Seems to be missing from canonical folder
        %         run_colorchange(myp,str,xyz, mycolors, actcolors);
        
        view(90,10); axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
    elseif strcmp(P,'cerebellum') || strcmp(P,'amygdala') || strcmp(P,'hypothalamus') ...
            || strcmp(P,'thalamus') || strcmp(P,'midbrain') || strcmp(P,'caudate') ...
            || strcmp(P,'globus pallidus') || strcmp(P,'putamen') || strcmp(P,'nucleus accumbens') ...
            || strcmp(P,'hippocampus')
        % this uses addbrain and works with any of its keywords
        p = [];
        
        myp = addbrain(P);p = [p myp];
        run_colorchange(myp,str,xyz, mycolors, actcolors);
        
        view(90,10); axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
        
    else
        error('Must input mat surf file or img file to convert to surf')
    end
    
end % if ishandle

lighting gouraud
lightRestoreSingle(gca);
material dull
axis off
set(gcf,'Color','w')
%scn_export_papersetup(400);  % this will mess up movies!!

if doverbose
    disp('Finished!')
    disp('___________________________________________')
end

end % main function






function actcolor = get_actcolors(datavaluesets, refZ, poscm, negcm)

actcolor = map_data_to_colormap(datavaluesets, poscm, negcm, refZ);

end




function run_colorchange(myp, str, xyz, mycolors, actcolors)

set(myp,'FaceAlpha',1);


for i = 1:length(myp)
    p = myp(i);
    
    % get original color
    origcolor = get(p,'FaceColor');
    
    % color change
    eval(str);
    
    % change non-active back to original color
    vdat = get(p,'FaceVertexCData');
    wh = find(all(vdat == .5, 2));
    vdat(wh,:) = repmat(origcolor,length(wh),1);
    set(p,'FaceVertexCData',vdat);
    
end
p = myp;
lighting gouraud; lightRestoreSingle(gca);

end





