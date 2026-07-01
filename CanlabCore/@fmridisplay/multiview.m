function obj = multiview(obj, montagetype, varargin)
% multiview  Compose a canonical multi-panel display onto this fmridisplay object.
%
% Lays out one of the montage/surface figure compositions defined in
% canlab_results_fmridisplay (e.g. 'compact', 'compact2', 'full', 'multirow',
% the HCP / freesurfer surface sets, the subcortical layouts, ...) onto THIS
% object, then renders any existing blob layers onto the new panels (via the
% fmridisplay.montage / surface pull-in). This is the fmridisplay method behind
% the default fmri_data.montage appearance; fmridisplay.montage(obj) routes here
% with 'compact'.
%
% :Usage:
% ::
%
%     obj = multiview(obj)                    % default 'compact' (sag + axial combo)
%     obj = multiview(obj, 'full')            % montages + surfaces
%     obj = multiview(obj, 'compact2')
%     obj = multiview(obj, montagetype, ...)  % any canlab_results_fmridisplay type
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object (handle). Composition is added to it in place.
%
%   **montagetype:**
%        A canlab_results_fmridisplay montage/surface layout keyword. Default
%        'compact'. See `help canlab_results_fmridisplay` for the full list
%        (compact, compact2, compact3, full, full2, multirow, allslices,
%        blobcenters/regioncenters, full hcp, hcp inflated, freesurfer inflated,
%        subcortex *, ...).
%
% :Optional Inputs:
%
%   Any further arguments are passed through to canlab_results_fmridisplay.
%
% :Outputs:
%
%   **obj:**
%        The same handle, with the composed montages/surfaces registered and
%        existing blob layers rendered onto them.
%
% :Examples:
% ::
%
%     o2 = fmridisplay;
%     o2 = addblobs(o2, region(my_stat_image));
%     o2 = multiview(o2, 'full');            % montages + surfaces, blobs pulled in
%
% :See also:
%   - canlab_results_fmridisplay, montage, surface, addblobs, fmridisplay
%
% ..
%    2026 visualization overhaul
% ..

if nargin < 2 || isempty(montagetype), montagetype = 'compact'; end

% Delegate the figure/axes composition to canlab_results_fmridisplay, passing
% THIS (handle) object so it composes onto it rather than building a new one
% (see its do_setup_display logic). 'noblobs' skips adding new blobs and
% 'noremove' keeps any EXISTING blob layers on the object (canlab_results_
% fmridisplay otherwise clears them before adding). The existing layers are
% rendered onto the new montages/surfaces by the fmridisplay.montage / surface
% pull-in as the panels are created. fmridisplay is a handle class, so obj is
% updated in place; we also capture the return for clarity.
obj = canlab_results_fmridisplay([], obj, montagetype, 'noblobs', 'noremove', varargin{:});

end
