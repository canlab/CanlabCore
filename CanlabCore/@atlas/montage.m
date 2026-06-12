function o2 = montage(obj, varargin)
% montage Display an atlas object on a standard slice montage.
%
% Convert the atlas object to a region object and forward to
% region/montage. By default, regions are colored with a standard color
% map; pass 'sourcespace' to override the source-space description used
% for surface projection.
%
% :Usage:
% ::
%
%     montage(obj, [optional inputs])
%
% Takes all optional inputs to canlab_results_fmridisplay and the
% addblobs method.
%
% :Inputs:
%
%   **obj:**
%        An atlas-class object.
%
% :Optional Inputs:
%
%   See help canlab_results_fmridisplay for the full list. Common
%   options include:
%
%   **o2:**
%        An existing fmridisplay object, with no keyword strings.
%
%   **'symmetric' [default]:**
%        Mirror left/right blobs with same colors. See
%        match_colors_left_right.
%
%   **'nosymmetric':**
%        Standard color map, no L/R color-matching for symmetry. Call
%        with 'nosymmetric' to match colors with wedge plots.
%
%   **'noblobs':**
%        Do not display blobs.
%
%   **'nooutline':**
%        Do not display blob outlines.
%
%   **'addmontages':**
%        When entering existing fmridisplay obj, add new montages.
%
%   **'noremove':**
%        Do not remove current blobs when adding new ones.
%
%   **'outlinecolor':**
%        Followed by new outline color.
%
%   **'splitcolor':**
%        Followed by a 4-cell new split colormap colors (help
%        fmridisplay or edit code for defaults as example).
%
%   **'montagetype':**
%        'full' for full montages of axial and saggital slices.
%        'full hcp' for full montage, but with surfaces and volumes from
%        HCP data.
%        'compact' [default] for single-figure parasagittal and axial slices.
%        'compact2': like 'compact', but fewer axial slices.
%        'multirow': followed by number of rows, e.g.,
%        o2 = canlab_results_fmridisplay([], 'multirow', 2).
%
%   **'noverbose':**
%        Suppress verbose output, good for scripts/publish to html, etc.
%
%   **'overlay':**
%        Specify anatomical image for montage (not surfaces), followed
%        by image name, e.g., o2 = canlab_results_fmridisplay([], ...
%        'overlay', 'icbm152_2009_symmetric_for_underlay.img').
%        The default brain for overlays is based on Keuken et al. 2014.
%        For legacy SPM8 single subject, enter as arguments: 'overlay',
%        which('SPM8_colin27T1_seg.img').
%
%   **'indexmap':**
%        Followed by a colormap (n x 3) to use as the indexed colormap
%        for atlas regions. Defaults to scn_standard_colors(num_regions).
%
%   **'sourcespace':**
%        Followed by a source-space identifier used by the surface
%        rendering routines (e.g., 'MNI152NLin2009cAsym'). Defaults to
%        obj.space_description.
%
% Other inputs to addblobs (fmridisplay method) are allowed, e.g.,
% 'cmaprange', [-2 2], 'trans'.
%
% :Outputs:
%
%   **o2:**
%        An fmridisplay object with the atlas regions rendered on slices
%        and (optionally) surfaces.
%
% :Examples:
% ::
%
%     o2 = canlab_results_fmridisplay([], 'noverbose');
%     o2 = montage(r, o2);      % symmetric colors left/right
%     o2 = removeblobs(o2);
%     o2 = montage(r, o2, 'map');
%
% :See also:
%   - region/montage
%   - canlab_results_fmridisplay
%   - fmridisplay
%   - atlas2region
%
% ..
%    Programmer Notes:
%    Tor Wager, Feb 2018. This function simply invokes the region montage
%    method.
% ..


cmap = [];
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'indexmap'
                cmap = varargin{i+1};
        end
    end
end

if isempty(cmap)
    nregions = num_regions(obj);
    cmap = scn_standard_colors(nregions);
    cmap = cell2mat(cmap');
    varargin = [varargin, 'indexmap', cmap];
end

wh = find(strcmp(varargin,'sourcespace'));
if ~isempty(wh)
    sourcespace = varargin{wh + 1};
    varargin{wh+1} = [];
    varargin{wh} = [];
else
    sourcespace = obj.space_description;
end

r = atlas2region(obj);

o2 = montage(r, 'sourcespace', sourcespace, varargin{:});


end % function


