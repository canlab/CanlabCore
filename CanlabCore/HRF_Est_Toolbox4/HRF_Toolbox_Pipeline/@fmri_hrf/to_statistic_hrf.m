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
%             statistic_hrf has t in .dat with .p computed under the
%             appropriate df.
%   'TStat'   already-computed t-statistic image. If provided, used directly
%             and 'SE' is ignored.
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
    t_obj = statistic_image(fmri_data(Hb));
    t_obj.dat = t_dat;
    t_obj.type = 'T';
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
