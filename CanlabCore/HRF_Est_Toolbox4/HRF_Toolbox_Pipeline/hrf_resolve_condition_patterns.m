function specs = hrf_resolve_condition_patterns(available_conditions, condition_spec, varargin)
%HRF_RESOLVE_CONDITION_PATTERNS Resolve exact/wildcard/regex conditions.
%
% specs = hrf_resolve_condition_patterns(available_conditions, condition_spec)
%
% condition_spec can be a numeric index, string, string array, or cellstr.
% Exact names are preferred. If no exact condition matches, wildcard
% patterns with * or ? are supported. Regex can be supplied as
% 'regex:<pattern>', '/<pattern>/', or a string containing common regex
% operators such as .*, |, [], (), ^, $, or +.
%
% The output is a struct array with indices and matched condition names. If
% one spec matches multiple conditions, downstream callers should average
% over specs(i).indices and display specs(i).display_label.

p = inputParser;
p.addRequired('available_conditions', @(x) iscell(x) || isstring(x));
p.addRequired('condition_spec', @(x) isempty(x) || isnumeric(x) || iscell(x) || isstring(x) || ischar(x));
p.addParameter('DefaultMode', 'each', @(x) ischar(x) || isstring(x));
p.addParameter('AllowNumeric', true, @(x) islogical(x) || isnumeric(x));
p.parse(available_conditions, condition_spec, varargin{:});
opts = p.Results;

available = cellstr(string(available_conditions));
default_mode = lower(char(opts.DefaultMode));

if isempty(condition_spec)
    specs = local_default_specs(available, default_mode);
    return
end

if isnumeric(condition_spec)
    if ~logical(opts.AllowNumeric)
        error('Numeric condition indices are not allowed here.');
    end
    specs = local_numeric_specs(available, condition_spec);
    return
end

tokens = cellstr(string(condition_spec));
tokens = tokens(~cellfun(@isempty, tokens));
if isempty(tokens)
    specs = local_default_specs(available, default_mode);
    return
end

specs = repmat(local_empty_spec(), 0, 1);
for i = 1:numel(tokens)
    specs(end + 1, 1) = local_token_spec(available, strtrim(tokens{i})); %#ok<AGROW>
end
end

function specs = local_default_specs(available, default_mode)
switch default_mode
    case 'each'
        specs = local_numeric_specs(available, 1:numel(available));
    case 'first'
        specs = local_numeric_specs(available, 1);
    case 'all'
        specs = local_group_spec(available, 1:numel(available), 'all', 'all');
    otherwise
        error('Unknown DefaultMode: %s. Use ''each'', ''first'', or ''all''.', default_mode);
end
end

function specs = local_numeric_specs(available, idx)
idx = idx(:)';
if any(idx < 1) || any(idx > numel(available)) || any(mod(idx, 1) ~= 0)
    error('Condition index out of range.');
end
specs = repmat(local_empty_spec(), numel(idx), 1);
for i = 1:numel(idx)
    specs(i) = local_group_spec(available, idx(i), available{idx(i)}, 'index');
end
end

function spec = local_token_spec(available, token)
idx = find(strcmp(available, token));
match_type = 'exact';

if isempty(idx)
    [pattern, match_type] = local_pattern(token);
    idx = find(~cellfun(@isempty, regexp(available, pattern, 'once')));
end

if isempty(idx)
    error('Condition pattern "%s" did not match any available condition. Available: %s', ...
        token, strjoin(available, ', '));
end

spec = local_group_spec(available, idx, token, match_type);
end

function [pattern, match_type] = local_pattern(token)
if startsWith(token, 'regex:')
    pattern = token(numel('regex:') + 1:end);
    match_type = 'regex';
elseif numel(token) >= 2 && startsWith(token, '/') && endsWith(token, '/')
    pattern = token(2:end - 1);
    match_type = 'regex';
elseif local_looks_like_regex(token)
    pattern = token;
    match_type = 'regex';
elseif contains(token, '*') || contains(token, '?')
    pattern = regexptranslate('wildcard', token);
    match_type = 'wildcard';
else
    pattern = ['^' regexptranslate('escape', token) '$'];
    match_type = 'exact';
end
end

function tf = local_looks_like_regex(token)
tf = contains(token, '.*') || contains(token, '|') || contains(token, '[') || ...
    contains(token, ']') || contains(token, '(') || contains(token, ')') || ...
    contains(token, '^') || contains(token, '$') || contains(token, '+');
end

function spec = local_group_spec(available, idx, label, match_type)
spec = local_empty_spec();
spec.pattern = char(label);
spec.match_type = char(match_type);
spec.indices = idx(:)';
spec.matched_conditions = available(spec.indices);
spec.label = char(label);
if isscalar(spec.indices)
    spec.display_label = spec.matched_conditions{1};
else
    spec.display_label = sprintf('%s (mean of %d: %s)', ...
        char(label), numel(spec.indices), strjoin(spec.matched_conditions, ', '));
end
end

function spec = local_empty_spec()
spec = struct('pattern', '', 'match_type', '', 'indices', [], ...
    'matched_conditions', {{}}, 'label', '', 'display_label', '');
end
