function [map, info] = assign_vals_to_atlas(atlas_obj, roi_names, vals, varargin)
% assign_vals_to_atlas  Project per-parcel values onto a brain map.
%
% Builds a brain image where every voxel in parcel i is assigned the value
% supplied for that parcel. Returns a statistic_image (default) so the
% standard .threshold / .montage / .table chain works, or an fmri_data.
%
% Generalizes the bespoke `assign_vals` helper in the Sun et al. workflow.
% Used internally by fmri_data.rsa_parcelwise.
%
% Usage
% -----
%   % Minimal: assign a t-value per parcel name
%   map = assign_vals_to_atlas(atlas, {'Ctx_V1_L','Ctx_V1_R'}, [3.1, 2.7]);
%
%   % All parcels in atlas order (roi_names = [] uses atlas.labels)
%   map = assign_vals_to_atlas(atlas, [], t_vals_per_parcel);
%
%   % Full statistic_image with p-values for thresholding
%   map = assign_vals_to_atlas(atlas, [], t_vals, ...
%       'p_vals', p_vals, 'output_type', 'statistic_image');
%   montage(threshold(map, 0.05, 'fdr'));   % function syntax: statistic_image
%                                           % has a `threshold` property AND method
%
% Inputs
% ------
%   atlas_obj  atlas object (the spatial template)
%   roi_names  cellstr of parcel labels to assign (must match atlas.labels).
%              Pass [] to assign all atlas parcels in label order (then
%              numel(vals) must equal num_regions(atlas)).
%   vals       numeric vector, one value per roi_name (or per parcel).
%
% Optional name-value
% -------------------
%   'output_type'  'statistic_image' (default) | 'fmri_data'
%   'p_vals'       numeric vector matching vals -- per-parcel p-values. When
%                  supplied (statistic_image output), fills .p so thresholding
%                  works. Default [] (sig set to all-true).
%   'fill'         value for voxels not in any assigned parcel. Default 0.
%   'dat_descrip'  string description stored on the output. Default ''.
%
% Outputs
% -------
%   map   statistic_image or fmri_data in atlas space
%   info  struct: .assigned_parcels, .n_assigned, .missing_names

p = inputParser;
p.addParameter('output_type', 'statistic_image', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'statistic_image','fmri_data'}));
p.addParameter('p_vals',      [],   @(x) isempty(x) || isnumeric(x));
p.addParameter('fill',        0,    @isnumeric);
p.addParameter('dat_descrip', '',   @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

if ~isa(atlas_obj, 'atlas')
    error('assign_vals_to_atlas:notAtlas', 'First argument must be an atlas object.');
end

labels = atlas_obj.labels;
n_parcels = num_regions(atlas_obj);

% Resolve roi_names -> parcel indices
if isempty(roi_names)
    if numel(vals) ~= n_parcels
        error('assign_vals_to_atlas:valCountMismatch', ...
            'roi_names empty implies all %d parcels, but numel(vals)=%d.', ...
            n_parcels, numel(vals));
    end
    parcel_idx = (1:n_parcels)';
    use_names  = labels(:);
else
    roi_names = cellstr(roi_names);
    if numel(vals) ~= numel(roi_names)
        error('assign_vals_to_atlas:valCountMismatch', ...
            'numel(roi_names)=%d but numel(vals)=%d.', numel(roi_names), numel(vals));
    end
    parcel_idx = nan(numel(roi_names), 1);
    for i = 1:numel(roi_names)
        hit = find(strcmp(labels, roi_names{i}), 1, 'first');
        if ~isempty(hit), parcel_idx(i) = hit; end
    end
    use_names = roi_names(:);
end

missing = isnan(parcel_idx);
info = struct();
info.missing_names   = use_names(missing);
if any(missing)
    warning('assign_vals_to_atlas:missingParcels', ...
        '%d roi_name(s) not found in atlas labels: %s', ...
        sum(missing), strjoin(info.missing_names, ', '));
end

% Build reference in atlas space; .dat holds integer parcel codes.
% Wrap in evalc to suppress the int32->single bit-rate conversion chatter.
[~, ref] = evalc('fmri_data(atlas_obj, ''noverbose'')');
codes = round(ref.dat(:, 1));      % integer parcel code per voxel

out_dat = opt.fill * ones(size(codes));
p_dat   = ones(size(codes));       % default p=1 (non-sig) where unassigned

have_p = ~isempty(opt.p_vals);
if have_p && numel(opt.p_vals) ~= numel(vals)
    error('assign_vals_to_atlas:pValCountMismatch', ...
        'numel(p_vals)=%d but numel(vals)=%d.', numel(opt.p_vals), numel(vals));
end

assigned = {};
for i = 1:numel(parcel_idx)
    if isnan(parcel_idx(i)), continue; end
    vmask = (codes == parcel_idx(i));
    if ~any(vmask), continue; end
    out_dat(vmask) = vals(i);
    if have_p, p_dat(vmask) = opt.p_vals(i); end
    assigned{end+1} = use_names{i}; %#ok<AGROW>
end
info.assigned_parcels = assigned;
info.n_assigned       = numel(assigned);

% Build the output object
if strcmpi(opt.output_type, 'fmri_data')
    map = ref;
    map.dat = out_dat;
    if ~isempty(char(opt.dat_descrip)), map.dat_descrip = char(opt.dat_descrip); end
else
    % statistic_image: start from the fmri_data reference and cast
    map = statistic_image('dat', out_dat, 'volInfo', ref.volInfo, ...
        'p', p_dat, 'type', 'generic', 'removed_voxels', ref.removed_voxels);
    map.dat = out_dat;
    map.p   = p_dat;
    % sig: if p_vals supplied, default sig at p<.05; else everything assigned is "sig"
    if have_p
        map.sig = p_dat < 0.05 & out_dat ~= opt.fill;
    else
        map.sig = (out_dat ~= opt.fill);
    end
    if ~isempty(char(opt.dat_descrip)), map.dat_descrip = char(opt.dat_descrip); end
end

end
