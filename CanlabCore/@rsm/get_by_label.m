function out = get_by_label(obj, query, varargin)
% get_by_label  Look up an rsm in an array by its source/parcel label.
%
% Searches the .source field of each rsm in the array for a match against
% the query string. Match modes (exact / contains / regex) selectable.
%
% Typical use case: a parcelwise compute_rsm call returns
%       R_per_parcel = [1 x nParcels] array of rsm objects
% each tagged with .source = 'parcel:<atlas_label>'. Looking up by index
% requires knowing the atlas label order; looking up by label does not.
%
% Usage
% -----
%   R_amy   = R_per_parcel.get_by_label('Amy')                  % contains (default)
%   R_acc   = R_per_parcel.get_by_label('CG_L', 'mode','exact')
%   R_sii   = R_per_parcel.get_by_label('SII\+?', 'mode','regex')
%
% If multiple parcels match, returns the full sub-array. If none match,
% errors with a list of available labels.
%
% Inputs
% ------
%   obj      rsm or array of rsm
%   query    char/string to match against .source
%   varargin:
%       'mode'         'contains' (default) | 'exact' | 'regex'
%       'ignore_case'  logical (default true)
%       'strip_prefix' logical (default true) -- ignore the 'parcel:' /
%                      'searchlight:' prefix when matching
%
% Output
% ------
%   out      Matching rsm (scalar if 1 match) or [1 x m] sub-array

p = inputParser;
p.addParameter('mode',          'contains', @(x) ischar(x) || isstring(x));
p.addParameter('ignore_case',   true,       @(x) islogical(x) || isnumeric(x));
p.addParameter('strip_prefix',  true,       @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

query = char(query);

n = builtin('numel', obj);
labels = cell(n, 1);
for i = 1:n
    s = char(obj(i).source);
    if opt.strip_prefix
        % Drop a leading 'kind:' prefix (parcel:, searchlight:, design:)
        c = strfind(s, ':');
        if ~isempty(c) && c(1) < numel(s)
            s = s(c(1)+1:end);
        end
    end
    labels{i} = s;
end

q = query;
if opt.ignore_case
    labels_cmp = lower(labels);
    q          = lower(q);
else
    labels_cmp = labels;
end

mode = lower(char(opt.mode));
switch mode
    case 'exact'
        hits = strcmp(labels_cmp, q);
    case 'contains'
        hits = contains(labels_cmp, q);
    case 'regex'
        hits = ~cellfun('isempty', regexp(labels_cmp, q, 'once'));
    otherwise
        error('rsm:get_by_label:badMode', ...
            'mode must be ''contains'', ''exact'', or ''regex''; got %s.', opt.mode);
end

idx = find(hits);
if isempty(idx)
    % Helpful error: list available labels
    preview = labels;
    if numel(preview) > 12, preview = [preview(1:10); {'...'}; preview(end-1:end)]; end
    error('rsm:get_by_label:noMatch', ...
        ['No rsm matched query "%s" (mode=%s). Available labels:\n  %s'], ...
        query, mode, strjoin(preview, ', '));
end

out = obj(idx);
% Return scalar when exactly one match (matches user expectation)
if numel(idx) == 1
    out = out(1);
end

end
