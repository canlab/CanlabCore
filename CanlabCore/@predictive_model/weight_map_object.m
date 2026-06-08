function [obj, si] = weight_map_object(obj, reference, varargin)
% weight_map_object  Build the model weight map and cache it on the model.
%
% Given a reference image object defined in the SAME voxel space the model
% was trained in, this maps the model's coefficient vector back into voxel
% space, wraps it as a @statistic_image, and STORES it on the model at
% pm.weights.weight_obj. After this, montage(pm) / surface(pm) /
% orthviews(weight_image(pm)) work with no further arguments.
%
% This is the piece that lets the numeric-input wrappers (xval_SVM,
% xval_SVR, xval_lasso_brain, ...) become montage/surface-ready: those
% functions are handed X/Y/id matrices and never see an image object, so
% they cannot build the weight map themselves. Train with the wrapper, then
% attach the map by handing this method any reference image in the same
% space (the fmri_data you pulled X out of is the obvious choice).
%
% If the model has been bootstrapped or feature-selected, the resulting
% weight object REFLECTS that: bootstrap z/p and the FDR-significant mask are
% attached to the @statistic_image (.p / .sig), and feature-selected models
% are expanded back to the full reference space (zeros for dropped voxels).
% You can therefore threshold/plot the cached object directly.
%
% :Usage:
% ::
%     pm = weight_map_object(pm, reference_fmri_data);
%     pm = weight_map_object(pm, reference_fmri_data, 'use', 'thresh_fdr');
%     [pm, si] = weight_map_object(pm, reference_fmri_data);
%
% :Inputs:
%
%   **obj:**
%        a fitted @predictive_model (after fit / crossval, and optionally
%        bootstrap / select_features).
%
%   **reference:**
%        an @fmri_data / @image_vector / @statistic_image whose
%        volInfo + removed_voxels map the weight vector back into voxel
%        space. Must be in the same space the model was trained in.
%
% :Optional Inputs (name/value, forwarded to weight_image):
%
%   **'use':**
%        which weight vector to map ('w' default, 'thresh_fdr',
%        'boot_w_mean'); see weight_image.
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with pm.weights.weight_obj populated.
%
%   **si:**
%        (optional) the @statistic_image weight map that was cached, in case
%        you want it directly.
%
% :Examples:
% ::
%     % --- numeric wrapper, then attach the map ---
%     dat = load_image_set('DPSP_hotwarm', 'noverbose');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = xval_SVM(X, Y, id);              % trained on matrices, no image
%     pm = weight_map_object(pm, dat);      % attach the weight map
%     montage(pm);                          % now works with no source
%     surface(pm);
%
%     % --- bootstrapped weights carry through to the object ---
%     pm = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
%     pm = weight_map_object(pm, dat);      % p / sig now attached
%     si = pm.weights.weight_obj;
%     si = threshold(si, .05, 'fdr');
%     montage(region(si));
%
% :See also:
%   weight_image, set_weight_obj, montage, surface, xval_SVM, xval_SVR

    if nargin < 2 || isempty(reference)
        error('predictive_model:weight_map_object:NoReference', ...
            ['Provide a reference image object (fmri_data / image_vector / ' ...
             'statistic_image) in the same space the model was trained in.']);
    end

    if ~isobject(reference) || ~isprop(reference, 'volInfo')
        error('predictive_model:weight_map_object:BadReference', ...
            ['reference must be an fmri_data / image_vector / statistic_image ' ...
             'with a volInfo property; got %s.'], class(reference));
    end

    % Build the statistic_image. weight_image already (a) expands
    % feature-selected weights back to full space and (b) attaches bootstrap
    % p / fdr_sig when they are present on the model, so the cached object
    % automatically reflects bootstrapping / feature selection.
    si = weight_image(obj, reference, varargin{:});

    % Cache it on the model in the canonical slot.
    obj = set_weight_obj(obj, si);

    % Provenance.
    has_boot = isstruct(obj.weights) && isfield(obj.weights, 'p') ...
        && ~isempty(obj.weights.p);
    if isprop(obj, 'history') && iscell(obj.history)
        obj.history{end+1} = sprintf( ...
            'weight_map_object: cached weight @statistic_image (%d voxels%s)', ...
            size(si.dat, 1), conditional_str(has_boot, ', with bootstrap p/sig', ''));
    end
end


function s = conditional_str(tf, a, b)
    if tf, s = a; else, s = b; end
end
