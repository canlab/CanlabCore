function [obj, si] = weight_map_object(obj, source, varargin)
% weight_map_object  Build the model weight map and cache it on the model.
%
% Given a reference image object defined in the SAME voxel space the model
% was trained in, this maps the model's coefficient vector back into voxel
% space, wraps it as a @statistic_image, and STORES it on the model at
% pm.weights.weight_obj. After this, montage(pm) / surface(pm) /
% orthviews(pm.weights.weight_obj) work with no further arguments.
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
% You can therefore threshold/plot the returned object directly.
%
% The second output is the @statistic_image itself, so this method also
% serves as the "just give me the weight image" builder — without caching,
% use [~, si] = weight_map_object(pm, source).
%
% :Usage:
% ::
%     pm        = weight_map_object(pm, reference_fmri_data);
%     pm        = weight_map_object(pm, reference_fmri_data, 'use', 'thresh_fdr');
%     [pm, si]  = weight_map_object(pm, reference_fmri_data);
%     [~,  si]  = weight_map_object(pm, reference_fmri_data);   % just the image
%
% :Inputs:
%
%   **obj:**
%        a fitted @predictive_model (after fit / crossval, and optionally
%        bootstrap / select_features).
%
%   **source:**
%        an @fmri_data / @image_vector / @statistic_image whose
%        volInfo + removed_voxels map the weight vector back into voxel
%        space. Must be in the same space the model was trained in.
%
% :Optional Inputs (name/value):
%
%   **'use':**
%        name of the weight vector to map. One of:
%          'w'           obj.weights.w (default; the full coefficient vector)
%          'thresh_fdr'  obj.weights.thresh_fdr (post-bootstrap)
%          'boot_w_mean' obj.weights.boot_w_mean
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with pm.weights.weight_obj populated.
%
%   **si:**
%        the @statistic_image weight map: the weight vector mapped into
%        source's voxel space (.dat), source's volInfo / removed_voxels, and
%        — when bootstrap stats are present — the bootstrap p-values in .p and
%        the FDR-significant mask in .sig, so you can threshold/plot directly.
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
%     % --- just the statistic_image, no caching ---
%     [~, si] = weight_map_object(pm, dat);
%
% :See also:
%   set_weight_obj, montage, surface, bootstrap, statistic_image,
%   xval_SVM, xval_SVR

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'use', 'w');
    parse(p, varargin{:});
    use = p.Results.use;

    if nargin < 2 || isempty(source)
        error('predictive_model:weight_map_object:NoSource', ...
            ['Provide a reference image object (fmri_data / image_vector / ' ...
             'statistic_image) in the same space the model was trained in.']);
    end

    if isempty(obj.weights) || ~isfield(obj.weights, use) || isempty(obj.weights.(use))
        error('predictive_model:weight_map_object:NoWeights', ...
            'pm.weights.%s is empty. Call fit() (and optionally bootstrap()) first.', use);
    end

    w = obj.weights.(use);
    if ~isnumeric(w) || ~isvector(w)
        error('predictive_model:weight_map_object:BadShape', ...
            'pm.weights.%s must be a numeric vector; got %s of size %s.', ...
            use, class(w), mat2str(size(w)));
    end
    w = w(:);

    if ~isobject(source) || ~isprop(source, 'volInfo')
        error('predictive_model:weight_map_object:BadSource', ...
            ['source must be an fmri_data / image_vector / statistic_image ' ...
             'with a volInfo property; got %s.'], class(source));
    end

    % If the model dropped features at fit time, expand back to the
    % full source feature space (zeros for dropped features). Check
    % the no-expansion case first — only fall through to the omitted-
    % features expansion if the lengths actually require it, since
    % many wrappers produce full-length weight vectors with an
    % all-false omitted_features.
    n_src = size(source.dat, 1);
    omitted = obj.omitted_features;
    if ~islogical(omitted) && isnumeric(omitted)
        omitted = logical(omitted);
    end

    if numel(w) == n_src
        full_w = w;
    elseif islogical(omitted) && numel(omitted) == n_src ...
            && numel(w) == n_src - sum(omitted)
        full_w = zeros(n_src, 1);
        full_w(~omitted) = w;
    else
        error('predictive_model:weight_map_object:LengthMismatch', ...
            ['Weights length (%d) does not fit source feature count (%d).\n' ...
             '  omitted_features length = %d, sum = %d.\n' ...
             'Was the model trained on a different fmri_data than `source`, or ' ...
             'has `source` been re-masked/remove_empty''d since fitting?'], ...
            numel(w), n_src, numel(omitted), sum(omitted ~= 0));
    end

    % Build the statistic_image.
    si = statistic_image();
    si.dat            = full_w;
    si.volInfo        = source.volInfo;
    if isprop(source, 'removed_voxels')
        si.removed_voxels = source.removed_voxels;
    end
    si.dat_descrip = sprintf('predictive_model weights (algorithm=%s, fit_type=%s, %s)', ...
        char(string(obj.algorithm)), char(string(obj.fit_type)), use);

    % Attach bootstrap stats if available, using the same fit-then-fall-back
    % logic as for full_w. This is what makes the cached object reflect
    % bootstrapping: .p / .sig let you threshold the result directly.
    has_boot = false;
    if isfield(obj.weights, 'p') && ~isempty(obj.weights.p) ...
            && numel(obj.weights.p) == numel(w)
        si.p = bridge_to_source_space(obj.weights.p, n_src, omitted, 1);
        has_boot = true;
    end
    if isfield(obj.weights, 'fdr_sig') && ~isempty(obj.weights.fdr_sig) ...
            && numel(obj.weights.fdr_sig) == numel(w)
        si.sig = logical(bridge_to_source_space(obj.weights.fdr_sig, n_src, omitted, 0));
    end

    % Cache it on the model in the canonical slot.
    obj = set_weight_obj(obj, si);

    % Provenance.
    if isprop(obj, 'history') && iscell(obj.history)
        boot_note = '';
        if has_boot, boot_note = ', with bootstrap p/sig'; end
        obj.history{end+1} = sprintf( ...
            'weight_map_object: cached weight @statistic_image (%d voxels%s)', ...
            size(si.dat, 1), boot_note);
    end
end


function out = bridge_to_source_space(v, n_src, omitted, fill_value)
% Helper: map a per-feature vector v from the model's feature space
% into the source's n_src-voxel space. Uses the same precedence as
% the main full_w branch: direct match wins, then omitted_features
% expansion, otherwise return v as-is (best-effort).
    v = v(:);
    if numel(v) == n_src
        out = v;
    elseif islogical(omitted) && numel(omitted) == n_src ...
            && numel(v) == n_src - sum(omitted)
        out = repmat(fill_value, n_src, 1);
        out(~omitted) = v;
    else
        out = v;
    end
end
