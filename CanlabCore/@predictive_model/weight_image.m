function si = weight_image(obj, source, varargin)
% weight_image  Return the model weights as a @statistic_image.
%
% :Usage:
% ::
%     si = weight_image(pm, source_fmri_data);
%     si = weight_image(pm, source_fmri_data, 'use', 'thresh_fdr');
%
% Maps the model's coefficient vector back into voxel space using
% source's volInfo. By default returns obj.weights.w (the full
% coefficient vector); pass 'use' to return a different weight vector
% (e.g. the FDR-thresholded version after bootstrap()).
%
% If pm.weights also has bootstrap stats (.z, .p, .fdr_thr,
% .fdr_sig), those are attached to the statistic_image so you can
% threshold/plot directly:
%
%     pm = bootstrap(pm, X, Y, 'nboot', 1000);
%     si = weight_image(pm, my_fmri_data);
%     si = threshold(si, .05, 'fdr');
%     montage(region(si));
%
% :Inputs:
%   obj       a fitted @predictive_model
%   source    an @fmri_data / @image_vector / @statistic_image to take
%             volInfo + removed_voxels from
%
% :Optional Inputs (name/value):
%   'use'     name of the weight vector to map. One of:
%               'w'           obj.weights.w (default)
%               'thresh_fdr'  obj.weights.thresh_fdr (post-bootstrap)
%               'boot_w_mean' obj.weights.boot_w_mean

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'use', 'w');
    parse(p, varargin{:});
    use = p.Results.use;

    if isempty(obj.weights) || ~isfield(obj.weights, use) || isempty(obj.weights.(use))
        error('predictive_model:weight_image:NoWeights', ...
            'pm.weights.%s is empty. Call fit() (and optionally bootstrap()) first.', use);
    end

    w = obj.weights.(use);
    if ~isnumeric(w) || ~isvector(w)
        error('predictive_model:weight_image:BadShape', ...
            'pm.weights.%s must be a numeric vector; got %s of size %s.', ...
            use, class(w), mat2str(size(w)));
    end
    w = w(:);

    if ~isobject(source) || ~isprop(source, 'volInfo')
        error('predictive_model:weight_image:BadSource', ...
            'source must be an fmri_data / image_vector / statistic_image with a volInfo property.');
    end

    % If the model dropped features at fit time, expand back to the
    % full source feature space (zeros for dropped features).
    n_src = size(source.dat, 1);
    if islogical(obj.omitted_features) && numel(obj.omitted_features) == n_src
        full_w = zeros(n_src, 1);
        full_w(~obj.omitted_features) = w;
    elseif numel(w) == n_src
        full_w = w;
    else
        error('predictive_model:weight_image:LengthMismatch', ...
            'weights length %d does not match source feature count %d (and omitted_features is the wrong shape to bridge).', ...
            numel(w), n_src);
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

    % Attach bootstrap stats if available.
    if isfield(obj.weights, 'p') && ~isempty(obj.weights.p) ...
            && numel(obj.weights.p) == numel(w)
        if islogical(obj.omitted_features) && numel(obj.omitted_features) == n_src
            full_p = ones(n_src, 1);
            full_p(~obj.omitted_features) = obj.weights.p;
        else
            full_p = obj.weights.p;
        end
        si.p = full_p;
    end
    if isfield(obj.weights, 'fdr_sig') && ~isempty(obj.weights.fdr_sig) ...
            && numel(obj.weights.fdr_sig) == numel(w)
        if islogical(obj.omitted_features) && numel(obj.omitted_features) == n_src
            full_sig = false(n_src, 1);
            full_sig(~obj.omitted_features) = obj.weights.fdr_sig;
        else
            full_sig = logical(obj.weights.fdr_sig);
        end
        si.sig = full_sig;
    end
end
