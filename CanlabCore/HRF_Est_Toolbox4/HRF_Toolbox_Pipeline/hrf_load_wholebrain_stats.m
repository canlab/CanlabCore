function wholebrain = hrf_load_wholebrain_stats(prefix, varargin)
%HRF_LOAD_WHOLEBRAIN_STATS Rebuild HRF statistic_image objects from files.
%
% wholebrain = hrf_load_wholebrain_stats(prefix)
%
% Reads files written by hrf_fit_wholebrain_stats:
%   <prefix>_beta.nii
%   <prefix>_t.nii
%   <prefix>_se.nii      optional for older outputs
%   <prefix>_p.nii       optional for older outputs
%   <prefix>_metadata.csv

p = inputParser;
p.addRequired('prefix', @(x) ischar(x) || isstring(x));
p.addParameter('NoVerbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(prefix, varargin{:});
opts = p.Results;

prefix = char(prefix);
paths = local_paths(prefix);
if exist(paths.beta, 'file') ~= 2
    error('Missing beta image: %s', paths.beta);
end
if exist(paths.t, 'file') ~= 2
    error('Missing T image: %s', paths.t);
end

load_args = {};
if logical(opts.NoVerbose)
    load_args = {'noverbose'};
end

b_obj = statistic_image(fmri_data(paths.beta, load_args{:}));
b_obj.type = 'FIR HRF beta';
t_obj = statistic_image(fmri_data(paths.t, load_args{:}));
t_obj.type = 'T';

metadata_table = local_read_metadata(paths.metadata);
[labels, N, dfe, TR] = local_metadata_fields(metadata_table, size(b_obj.dat, 2));

if ~isempty(labels)
    b_obj.image_labels = labels;
    t_obj.image_labels = labels;
end

if ~isempty(N)
    b_obj.N = N;
    t_obj.N = N;
end
if ~isempty(dfe)
    b_obj.dfe = dfe;
    t_obj.dfe = dfe;
end

if exist(paths.se, 'file') == 2
    se_obj = fmri_data(paths.se, load_args{:});
    b_obj.ste = se_obj.dat;
    t_obj.ste = se_obj.dat;
elseif ~isempty(b_obj.dat) && ~isempty(t_obj.dat)
    ste = b_obj.dat ./ t_obj.dat;
    ste(~isfinite(ste)) = NaN;
    b_obj.ste = ste;
    t_obj.ste = ste;
end

if exist(paths.p, 'file') == 2
    p_obj = fmri_data(paths.p, load_args{:});
    b_obj.p = p_obj.dat;
    t_obj.p = p_obj.dat;
elseif ~isempty(dfe)
    t_obj.dat = double(t_obj.dat);
    p_dat = 2 * (1 - tcdf(abs(t_obj.dat), dfe));
    p_dat(p_dat == 0) = eps;
    b_obj.p = p_dat;
    t_obj.p = p_dat;
end

if isempty(b_obj.p_type), b_obj.p_type = 'Two-tailed P-values from HRF T map'; end
if isempty(t_obj.p_type), t_obj.p_type = 'Two-tailed P-values from HRF T map'; end
if isempty(b_obj.sig), b_obj.sig = true(size(b_obj.dat)); end
if isempty(t_obj.sig), t_obj.sig = true(size(t_obj.dat)); end

wholebrain = struct();
wholebrain.b = b_obj;
wholebrain.t = t_obj;
wholebrain.metadata_table = metadata_table;
if ~isempty(metadata_table) && any(strcmp('condition', metadata_table.Properties.VariableNames))
    wholebrain.conditions = unique(cellstr(string(metadata_table.condition)), 'stable');
else
    wholebrain.conditions = {};
end
wholebrain.TR = TR;
wholebrain.dfe = dfe;
wholebrain.N = N;
wholebrain.paths = local_existing_paths(paths);
wholebrain.reused = true;
end

function paths = local_paths(prefix)
paths = struct();
paths.beta = [prefix '_beta.nii'];
paths.t = [prefix '_t.nii'];
paths.se = [prefix '_se.nii'];
paths.p = [prefix '_p.nii'];
paths.metadata = [prefix '_metadata.csv'];
paths.t_thresholded = [prefix '_t_thresh.nii'];
end

function metadata_table = local_read_metadata(metadata_file)
if exist(metadata_file, 'file') == 2
    metadata_table = readtable(metadata_file, 'TextType', 'string');
else
    metadata_table = table();
end
end

function [labels, N, dfe, TR] = local_metadata_fields(metadata_table, n_images)
labels = {};
N = [];
dfe = [];
TR = [];
if isempty(metadata_table)
    return
end

if any(strcmp('image_label', metadata_table.Properties.VariableNames)) && height(metadata_table) == n_images
    labels = cellstr(string(metadata_table.image_label));
end
if any(strcmp('N', metadata_table.Properties.VariableNames))
    N = metadata_table.N(1);
end
if any(strcmp('dfe', metadata_table.Properties.VariableNames))
    dfe = metadata_table.dfe(1);
end
if any(strcmp('TR', metadata_table.Properties.VariableNames))
    TR = metadata_table.TR(1);
end
end

function paths = local_existing_paths(paths)
names = fieldnames(paths);
for i = 1:numel(names)
    if exist(paths.(names{i}), 'file') ~= 2
        paths = rmfield(paths, names{i});
    end
end
end
