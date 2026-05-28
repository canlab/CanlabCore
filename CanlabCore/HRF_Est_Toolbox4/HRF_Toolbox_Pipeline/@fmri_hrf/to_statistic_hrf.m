function Ht = to_statistic_hrf(Hb, varargin)
% to_statistic_hrf - Build a paired statistic_hrf from an fmri_hrf.
%
% Usage
% -----
%   Ht = to_statistic_hrf(Hb, 'SE', se_obj)
%   Ht = to_statistic_hrf(Hb, 'TStat', t_obj)
%
% Inputs
% ------
%   Hb        fmri_hrf, holding beta values.
%   'SE'      statistic_image (or fmri_data) holding standard errors aligned
%             with Hb. If provided, t = beta ./ se is computed; the resulting
%             statistic_hrf carries that t in .dat.
%   'TStat'   already-computed t-statistic image. If provided, used directly.
%
% Returns a statistic_hrf carrying the same HRF metadata as Hb.

p = inputParser;
p.addParameter('SE', [], @(x) isempty(x) || isa(x, 'image_vector'));
p.addParameter('TStat', [], @(x) isempty(x) || isa(x, 'image_vector'));
p.parse(varargin{:});
opts = p.Results;

if isempty(opts.TStat) && isempty(opts.SE)
    error('fmri_hrf:to_statistic_hrf:NoStat', ...
        'Provide either an SE image (to derive t from beta/SE) or a TStat image.');
end

% Build a statistic_image with the t-value data, lifting the image_vector
% layout (voxel space, removed_voxels, etc.) from Hb. We do this by hand
% rather than via statistic_image(fmri_data(Hb)), because fmri_data's
% constructor expects an image filename as its first arg and would
% misinterpret an existing object.
if ~isempty(opts.TStat)
    t_obj = opts.TStat;
else
    se = opts.SE;
    if ~isequal(size(Hb.dat), size(se.dat))
        error('fmri_hrf:to_statistic_hrf:ShapeMismatch', ...
            'SE shape %s does not match beta shape %s.', ...
            mat2str(size(se.dat)), mat2str(size(Hb.dat)));
    end
    t_dat = Hb.dat ./ se.dat;
    t_obj = local_build_t_image(Hb, t_dat);
end

Ht = statistic_hrf(t_obj, ...
    'MetadataTable', Hb.metadata_table, ...
    'Subject', Hb.subject, ...
    'RunLabel', Hb.run_label, ...
    'ModelName', Hb.model_name, ...
    'Conditions', Hb.conditions, ...
    'TR', Hb.TR, ...
    'DesignMatrix', Hb.design_matrix, ...
    'DesignInfo', Hb.design_info, ...
    'SourcePaths', Hb.source_paths);
end


function si = local_build_t_image(Hb, t_dat)
% Construct a statistic_image with the t-values, copying the voxel-space
% layout from Hb (an fmri_hrf, which inherits its image_vector fields).
si = statistic_image();
shared = {'dat', 'volInfo', 'removed_voxels', 'removed_images', ...
    'space_defining_image_name', 'fullpath', 'files_exist', 'history', ...
    'image_names', 'source_notes'};
for i = 1:numel(shared)
    f = shared{i};
    if isprop(si, f) && isprop(Hb, f)
        try
            si.(f) = Hb.(f);
        catch
        end
    end
end
si.dat = t_dat;
si.type = 'T';
if isprop(si, 'sig') && isempty(si.sig)
    si.sig = true(size(si.dat));
end
end
