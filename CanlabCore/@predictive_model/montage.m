function varargout = montage(obj, varargin)
% montage  Montage of a predictive model's weight map.
%
% Thin delegate: build the weight @statistic_image (via weight_image) and
% call its montage. With 'regions', route through region() and montage the
% region object instead, which outlines contiguous significant clusters.
%
% :Usage:
% ::
%     montage(pm, source_fmri_data);
%     montage(pm, source_fmri_data, 'regions');
%     montage(pm, source_fmri_data, 'use', 'thresh_fdr');
%     montage(pm);                 % uses pm.weights.weight_obj if present
%
% :Inputs:
%
%   **obj:**
%        a fitted @predictive_model.
%
%   **source:**
%        (optional) an @fmri_data / @image_vector / @statistic_image whose
%        volInfo + removed_voxels map the weight vector back into voxel
%        space. May be omitted if obj.weights.weight_obj is already
%        populated (e.g. by a brain wrapper or fmri_data.predict 'newapi').
%
% :Optional Inputs (name/value, forwarded):
%
%   **'use':**
%        which weight vector to map ('w' default, 'thresh_fdr',
%        'boot_w_mean'); see weight_image.
%
%   **'regions':**
%        montage region(weight_image(...)) instead of the statistic_image
%        directly. Most useful after bootstrap() has set FDR significance.
%
%   Any other name/value pairs are passed through to the underlying
%   image_vector.montage / region.montage.
%
% :Outputs:
%
%   **fig_handle:**
%        the figure handle returned by the underlying montage call.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = crossval(pm, X, Y);
%     pm = bootstrap(pm, X, Y, 'nboot', 1000);
%     montage(pm, dat);                      % weight map on slices
%     montage(pm, dat, 'use', 'thresh_fdr', 'regions');  % FDR clusters
%
% :See also:
%   weight_image, weight_map_object, surface, region, image_vector.montage

    [si, passthrough, do_regions] = predictive_model.weight_image_for_display(obj, varargin{:});

    if do_regions
        r = region(si);
        if nargout > 0
            [varargout{1:nargout}] = montage(r, passthrough{:});
        else
            montage(r, passthrough{:});
        end
    else
        if nargout > 0
            [varargout{1:nargout}] = montage(si, passthrough{:});
        else
            montage(si, passthrough{:});
        end
    end
end
