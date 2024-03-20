function o2 = montage(obj, varargin)
% This function displays a region object on a standard slice montage
%
% :Usage:
% ::
%
%    montage(obj, [optional inputs])
%
% - takes all optional inputs to canlab_results_fmridisplay and addblobs method
%
% :Input:
%
%   **obj:**
%        a region object
%
% :Optional Inputs:  - see help canlab_results_fmridisplay
%
%   **o2***
%        An existing fmridisplay object, with no keyword strings
%
%   **'color'**
%        Followed by a single [r g b] triplet to define color for all regions. Passed
%        through to addblobs.
%
%   **'colormap'**
%        Map color according to statistic values/data entered in .Z field
%
%   **'colors'**
%        Followed by a cell array of [r g b] triplets, one color for each
%        region. Note: If 'symmetric' option is on, will reorganize color
%        to region assignments.
%        Note: if an n x 3 matrix is supplied instead of a cell array a
%        warning is thrown and automatic conversion to the correct format
%        is attempted.
%
%   **'indexmap'**
%        Interprets data as an indexmap and plots as a single 'blob'
%        mapping colors to data values.
%        Note: {'colors', cmap} and {'indexmap', cmap} may produce similar
%        results in certain cases, but {'colors', cmap} may result in some
%        overlap between adjacent regions due to an interpolation artifact 
%        while {'indexmap',cmap} doesn't suffer from this problem. This
%        makes the indexmap better for displaying parcellations on a brain.
%
%   **symmetric** 
%       Mirror left/right blobs with same colors
%       See match_colors_left_right
%
%   **nosymmetric** [default]
%       Standard color map, no L/R color-matching for symmetry 
%
%   **'noblobs':**
%        do not display blobs
%
%   **'nooutline':**
%        do not display blob outlines
%
%   **'nofigure':**
%        Create slices on existing figure; do not clear 'fmridisplay'-tagged figure
%
%   **'addmontages':**
%        when entering existing fmridisplay obj, add new montages
%
%   **'noremove':**
%        do not remove current blobs when adding new ones
%
%   **'outlinecolor:**
%        followed by new outline color
%
%   **'splitcolor':**
%        followed by 4-cell new split colormap colors (help fmridisplay or edit code for defaults as example)
%
%   **'montagetype':**
%        'full' for full montages of axial and sagg slices.
%
%        'full hcp' for full montage, but with surfaces and volumes from
%        HCP data
%
%        'full hcp inflated' for full montage with inflated surfaces from
%        HCP data
%
%        'compact' [default] for single-figure parasagittal and axials slices.
%
%        'compact2': like 'compact', but fewer axial slices.
%
%        'multirow': followed by number of rows
%           e.g., o2 = canlab_results_fmridisplay([], 'multirow', 2);
%
%        'regioncenters': Show one slice per region, centered on region center
%               Note: 'nozoom' option omits zooming-in on regions
%
%   **'noverbose':**
%        suppress verbose output, good for scripts/publish to html, etc.
%
%   **'overlay':**
%        specify anatomical image for montage (not surfaces), followed by
%        image name
%        e.g., o2 = canlab_results_fmridisplay([], 'overlay', 'icbm152_2009_symmetric_for_underlay.img')';
%
%         The default brain for overlays is based on Keuken et al. 2014
%         For legacy SPM8 single subject, enter as arguments:
%         'overlay', which('SPM8_colin27T1_seg.img')
%
% Other inputs to addblobs (fmridisplay method) are allowed, e.g., 'cmaprange', [-2 2], 'trans'
%
% See help fmridisplay
% e.g., 'color', [1 0 0]
%
%
% :Examples:
% -------------------------------------------------------------------------
% ::
% Example 1: 
% % Complete group analysis of a standard dataset
% % Do analysis and prep results region object:
%
%   img_obj = load_image_set('emotionreg');         % Load a dataset
%   t = ttest(img_obj, .005, 'unc');                % Do a group t-test
%   t = threshold(t, .005, 'unc', 'k', 10);         % Re-threshold with extent threshold of 10 contiguous voxels
%   r = region(t);                                  % Turn t-map into a region object with one element per contig region
%
%   Label regions and print a table:
%   [r, region_table, table_legend_text] =
%   autolabel_regions_using_atlas(r);               % Label regions. Can be skipped because 'table' below attempts to do this automatically
%   table(r);                                       % Print a table of results using new region names
%
%   Display montages in several styles:
%
%   montage(r)                                          % A unique color per blob (the default)
%   montage(r, 'colormap')                              % A color-map voxels according to statistic values/intensity
%   montage(r, 'color', [1 0 0])                        % All blobs in red
%   montage(r, 'colormap', 'full')                      % montage type is 'full', with surfaces
%   montage(r, 'colormap', 'regioncenters')             % montage type is 'regioncenters', with slices located at the center of each blob (zooms in on blobs)
%   
%   Use colormap, but with blue colors. see addblobs for more options:
%    montage(r, 'colormap', 'maxcolor', [0 0 1], 'mincolor', [.6 .5 .8], 'regioncenters');   
%
%  Pass out "o2", an fmridisplay object that has blobs registered, so you
%  can remove them and reuse them:
%
%  o2 = montage(r, 'regioncenters', 'colormap');
%  o2 = removeblobs(o2);
%  o2 = montage(r, o2, 'regioncenters', 'color', [1 0 0]);
%
% Example 2: 
% % Extend previous results by creating custom fmridisplay object
% % and adding blobs to that.
%
%   o2 = canlab_results_fmridisplay([], 'noverbose');
%   o2 = montage(r, o2);      % symmetric colors left/right
%   o2 = removeblobs(o2);
%   o2 = montage(r, o2, 'map');
% 
% Example 3: 
% % create a custom fmridisplay object and display regions from a standard
% parcellation (Glasser 2016 Nature cortical parcellation)
% 
%   create_figure('slices'); axis off
%   o2 = canlab_results_fmridisplay([], 'multirow', 2);
%   brighten(.6)
%   hcp152t1 = which('HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii');
%   r = region(fmri_data(hcp152t1), 'unique_mask_values');
%   o2 = montage(r, o2, 'wh_montages', 3:4);
%   o2 = montage(r(1:20), o2, 'wh_montages', 1:2, 'color', [1 .5 0]);
%
%   % Use a different montage type in canlab_results_fmridisplay:
%   o2 = montage(r, 'compact2', 'nosymmetric');
%
% % Plot one region blob per slice on a series of montages.  
%   o2 = montage(r, 'regioncenters', 'nosymmetric');

% edited: Tor, 1/2018, unique color option/default
%         Tor  7/2018 .added documentation

% Defaults and inputs
% -----------------------------------------------------------------------

methodtype = 'map'; % 'symmetric';
colortype = 'unique';
one_blob_per_slice = false; 
dozoom = true;
colors = scn_standard_colors(length(obj));
dofigure = true;
doredefinecolors = true;

if any(strcmp(varargin, 'map')), methodtype = 'map'; end
if any(strcmp(varargin, 'nosymmetric')), methodtype = 'map'; end
if any(strcmp(varargin, 'symmetric')), methodtype = 'symmetric'; end
if any(strcmp(varargin, 'old')), methodtype = 'old'; end
if any(strcmp(varargin, 'nofigure')), dofigure = false; end
if any(strcmp(varargin, 'indexmap'))
    colors = varargin{find(strcmp(varargin,'indexmap'))+1};
    colors = mat2cell(colors,repmat(1,1,length(colors)),3)';
    colortype = 'indexmap';
    
    if any(strcmp(varargin,'symmetric'))
        warning('''symmetric'' argument is ignored when providing an explicit indexed colormapping.');
    end
end

if any(strcmp(varargin, 'color'))
    colortype = 'solid'; 
    % No need to pass in color here, because varargin{:} is passed forward,
    % and interpreted by addblobs   
end

if any(strcmp(varargin, 'colormap')) 
    % Without any other arguments, this uses a colormap by default in
    % subfunctions, since no 'color' arg is passed forward
    
    colortype = 'solid'; 

end

if any(strcmp(varargin, 'colors'))
    % If we have passed in specific colors, don't redefine them
    wh = find(strcmp(varargin, 'colors'));
    colors = varargin{wh + 1}; varargin{wh} = []; varargin{wh + 1} = []; 
    if ~iscell(colors) && size(colors,2) == 3
        warning('''colors'' supplied as matrix, assuming this is a colormap n x 3 matrix, and converting accordingly');
        colors = mat2cell(colors,repmat(1,1,length(colors)),3)';
    elseif ~iscell(colors)
        error('colors must be a cell array of 3-vectors');
    end
    
    doredefinecolors = false;
end

if any(strcmp(varargin, 'regioncenters')), one_blob_per_slice = true; end
if any(strcmp(varargin, 'nozoom')), dozoom = false; end

% Initialize display if needed
% -----------------------------------------------------------------------

funhan = @(x) isa(x, 'fmridisplay');
whfmridisplay = cellfun(funhan, varargin);

o2 = [];

if any(whfmridisplay)
    wh = find(whfmridisplay);
    o2 = varargin{wh(1)};
end

if ~exist('o2', 'var') || ~isa(o2, 'fmridisplay')
    
    if dofigure && ~(one_blob_per_slice) % regioncenters will create a new figure anyway
        create_figure('fmridisplay'); axis off
    end
    % TW: I think the above is unnecessary because canlab_results_fmridisplay
    % will create a figure anyway if needed. and w/o this can pass in
    % nofigure control string. but needs refactoring to work better.
    
    % some inputs control which slices are shown: standard slices, or specific locations for these blobs
    % regioncenters is passed through, so it determines the type of montage to create:
    o2 = canlab_results_fmridisplay(obj, 'noverbose', varargin{:}, 'noblobs', 'nooutline');

end


% Define colors for 'unique' option
% And render blobs
% -------------------------------------------------------------------------
switch colortype
    
    case 'unique' % different color for each blob (some reuse)
        indx = find(strcmp(varargin,'indexmap'));
        if ~isempty(indx)
            varargin(indx:indx+1) = [];
        end
        
        switch methodtype
            case 'symmetric'
                colors = match_colors_left_right(obj, @(x1)(colors(1:x1)));
        end
        
        if one_blob_per_slice
            % One blob on each slice
            
            for i = 1:length(obj)
                
                o2 = addblobs(o2, obj(i), 'color', colors{i}, 'noverbose', varargin{:}, 'wh_montages', i);
               
                o2 = title_montage(o2, i, obj(i).shorttitle);

            end
            
            if dozoom
                wh = find(cellfun(@ischar,varargin));
                wh = wh(ismember(varargin(wh),{'axial','saggital','coronal'}));
                if any(wh)
                    zoom_in_on_regions(o2, obj, varargin{wh});
                else
                    zoom_in_on_regions(o2, obj, 'axial'); 
                end
            end
            
            drawnow
            
        else
            % All visible blobs on each slice
            
            if doredefinecolors
                % Redefine blob groups based on color categories, to reduce number of regions, combining those with same color too many function calls are bad news!! overloads handle graphics system?
                [region_groups, colors] = redefine_colors_based_on_groups(obj, colors);
                
                for i = 1:length(region_groups)
                    
                    o2 = addblobs(o2, region_groups{i}, 'color', colors{i}, 'noverbose', varargin{:});
                    drawnow
                    
                end
                
            else
                % no redefine
                
                for i = 1:length(obj)
                    
                    o2 = addblobs(o2, obj(i), 'color', colors{i}, 'noverbose', varargin{:});
                    drawnow
                    
                end
                
            end
            
        end
        
    case {'solid', 'indexmap'} % solid color, user-entered, or use colormap
        
        if one_blob_per_slice
            
            for i = 1:length(obj)
                
                o2 = addblobs(o2, obj(i), varargin{:}, 'wh_montages', i, 'noverbose');
                o2 = title_montage(o2, i, obj(i).shorttitle);
                
            end
            
            if dozoom
                zoom_in_on_regions(o2, obj, 'axial');  % hard-coded for now, could change orientation...must make flexible in canlab_results_fmridisplay
            end
            
            drawnow
            
        else
            
            % Just render blobs
            o2 = addblobs(o2, obj, varargin{:});
        end
        
    case 'old'
        
        montage_clusters([], obj, varargin{:})
        
end


end % function


function [region_groups, color_groups] = redefine_colors_based_on_groups(r, colors)

% regroup to reduce number of regions, combining those with same color
% Why? too many function calls are bad news!! overloads handle graphics system.

if iscell(colors)
    catcolors = cat(1, colors{:});
    
elseif size(colors, 2) ~= 3
    % Maybe wrong orientation
    catcolors = colors';
    disp('Bad input: colors is not cell. check code.')
    keyboard
    colors = {colors'};     % this will not work
else
    catcolors = colors;
    colors = {colors};
end

uniquecolors = unique(catcolors, 'rows');
ncolors = size(uniquecolors, 1);

region_groups = {};
color_groups = {};

for i = 1:ncolors
    
    mycolor = uniquecolors(i, :);
    cmatchfun = @(x) ~any(x - mycolor);
    
    wh_regions = cellfun(cmatchfun, colors);
    
    region_groups{i} = r(wh_regions);
    
    color_groups{i} = mycolor;
    
    % o2 = addblobs(o2, r(wh_regions), 'maxcolor', mycolor, 'mincolor', mycolor, 'wh_montages', 1:2);
    
end

end
