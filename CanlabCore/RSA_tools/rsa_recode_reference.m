function out = rsa_recode_reference(values, reference_levels, varargin)
% rsa_recode_reference  Collapse a metadata column to {reference, other} levels.
%
% For "shared-anchor" designs where every subject shares one or more
% reference level(s) of a factor but each has idiosyncratic other levels
% (e.g. all subjects have a "Left Face" bodysite plus one different "other"
% bodysite each), this recodes the column so the reference level(s) are kept
% verbatim and every other value collapses to a single label. The recoded
% column is consistent across subjects, so compute_rsm builds the same k
% conditions for everyone and cross-subject reliability / contrasts become
% well-defined.
%
% :Usage:
% ::
%     out = rsa_recode_reference(values, reference_levels, ['other_label', label])
%
% :Inputs:
%   **values:**
%        a metadata column (cellstr / string / categorical / numeric).
%
%   **reference_levels:**
%        the value(s) to keep verbatim (char, string, cellstr, or numeric).
%        Everything else collapses to the other_label.
%
% :Optional Inputs:
%   **'other_label':**
%        label for the collapsed non-reference values. Default 'Other'.
%
%   **'ignore_case':**
%        case-insensitive matching of reference levels. Default true.
%
% :Outputs:
%   **out:**
%        a cellstr column the same height as values, with reference levels
%        preserved and all others set to other_label.
%
% :Examples:
% ::
%     % DistractMap / AcceptMap: every subject has Left Face + one other site
%     T = distractmap_run.metadata_table;
%     T.bodysite_type = rsa_recode_reference(T.bodySite, 'Left Face', ...
%                                            'other_label', 'Other Body Site');
%     distractmap_run.metadata_table = T;
%
%     R = compute_rsm(distractmap_run, 'group_by', {'condition','bodysite_type'}, ...
%                     'subject_var','subject_id', 'level','session', ...
%                     'session_var','session_number', 'metric','spearman');
%     out = R.reliability_per_condition;
%
% :See also: compute_rsm, rsm.reliability_per_condition
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

% ---------- INPUT PARSER ----------
p = inputParser;
p.addRequired('values');
p.addRequired('reference_levels');
p.addParameter('other_label', 'Other', @(x) ischar(x) || isstring(x));
p.addParameter('ignore_case', true, @(x) islogical(x) || isnumeric(x));
p.parse(values, reference_levels, varargin{:});
other_label = char(p.Results.other_label);
ignore_case = logical(p.Results.ignore_case);

% ---------- Normalize inputs to cellstr ----------
vals_str = local_to_cellstr(values);
refs     = local_to_cellstr(reference_levels);

% ---------- Match ----------
if ignore_case
    is_ref = ismember(lower(vals_str), lower(refs));
else
    is_ref = ismember(vals_str, refs);
end

out = repmat({other_label}, numel(vals_str), 1);
out(is_ref) = vals_str(is_ref);

% Warn if a reference level never matched (likely a typo / case issue)
matched_refs = unique(vals_str(is_ref));
for i = 1:numel(refs)
    if ignore_case
        hit = any(strcmpi(matched_refs, refs{i}));
    else
        hit = any(strcmp(matched_refs, refs{i}));
    end
    if ~hit
        warning('rsa_recode_reference:refNotFound', ...
            'Reference level "%s" did not match any value in the column.', refs{i});
    end
end

end


function c = local_to_cellstr(x)
if iscell(x)
    c = cellfun(@(v) char(string(v)), x, 'UniformOutput', false);
    c = c(:);
elseif iscategorical(x)
    c = cellstr(x); c = c(:);
elseif isstring(x)
    c = cellstr(x); c = c(:);
elseif ischar(x)
    c = {x};
elseif isnumeric(x) || islogical(x)
    c = cellfun(@(v) char(string(v)), num2cell(x(:)), 'UniformOutput', false);
else
    error('rsa_recode_reference:badType', 'Unsupported input type: %s', class(x));
end
end
