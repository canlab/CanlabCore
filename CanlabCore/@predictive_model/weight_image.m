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
%
% :Outputs:
%
%   **si:**
%        a @statistic_image with the weight vector mapped into source's
%        voxel space (.dat), source's volInfo / removed_voxels, and — when
%        bootstrap stats are present — the bootstrap p-values in .p and the
%        FDR-significant mask in .sig, so you can threshold/plot directly.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = crossval(pm, X, Y);
%     pm = bootstrap(pm, X, Y, 'nboot', 1000);
%     si = weight_image(pm, dat);        % statistic_image with .p / .sig
%     si = threshold(si, .05, 'fdr');
%     montage(region(si));
%
% :See also:
%   montage, surface, bootstrap, statistic_image

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
        error('predictive_model:weight_image:LengthMismatch', ...
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
    % logic as for full_w.
    if isfield(obj.weights, 'p') && ~isempty(obj.weights.p) ...
            && numel(obj.weights.p) == numel(w)
        si.p = bridge_to_source_space(obj.weights.p, n_src, omitted, 1);
    end
    if isfield(obj.weights, 'fdr_sig') && ~isempty(obj.weights.fdr_sig) ...
            && numel(obj.weights.fdr_sig) == numel(w)
        si.sig = logical(bridge_to_source_space(obj.weights.fdr_sig, n_src, omitted, 0));
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
