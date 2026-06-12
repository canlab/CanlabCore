function varargout = surface(obj, varargin)
% surface  Cortical-surface rendering of a predictive model's weight map.
%
% Thin delegate: build the weight @statistic_image (via weight_map_object)
% and call its surface. With 'regions', route through region() and surface
% the region object instead.
%
% :Usage:
% ::
%     surface(pm, source_fmri_data);
%     surface(pm, source_fmri_data, 'use', 'thresh_fdr');
%     surface(pm);                 % uses pm.weights.weight_obj if present
%
% :Inputs:
%
%   **obj:**
%        a fitted @predictive_model.
%
%   **source:**
%        (optional) an @fmri_data / @image_vector / @statistic_image whose
%        volInfo + removed_voxels map the weight vector back into voxel
%        space. May be omitted if obj.weights.weight_obj is populated.
%
% :Optional Inputs (name/value, forwarded):
%
%   **'use':**
%        which weight vector to map ('w' default, 'thresh_fdr', ...).
%
%   **'regions':**
%        surface region(weight_map_object(...)) instead of the statistic_image.
%
%   Any other name/value pairs are passed through to the underlying
%   image_vector.surface / region.surface.
%
% :Outputs:
%
%   **all_surf_handles:**
%        surface handles returned by the underlying surface call (plus the
%        positive/negative cluster objects for the region path).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = crossval(pm, X, Y);
%     pm = bootstrap(pm, X, Y, 'nboot', 1000);
%     surface(pm, dat, 'use', 'thresh_fdr');
%
% :See also:
%   weight_map_object, montage, region, image_vector.surface

    [si, passthrough, do_regions] = predictive_model.weight_image_for_display(obj, varargin{:});

    if do_regions
        r = region(si);
        if nargout > 0
            [varargout{1:nargout}] = surface(r, passthrough{:});
        else
            surface(r, passthrough{:});
        end
    else
        if nargout > 0
            [varargout{1:nargout}] = surface(si, passthrough{:});
        else
            surface(si, passthrough{:});
        end
    end
end
