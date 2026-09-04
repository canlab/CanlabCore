function col = detect_subject_column(replicate_table)
% detect_subject_column  Find the subject-id column in an rsm.replicate_table.
%
% Returns the first matching column name (case-insensitive) from this priority
% list, or '' if none match:
%   subject_id, sub, subject, subjectid, subj, sid
%
% Used by the reliability methods so per-subject ICC is the default whenever
% a subject axis is present in the replicate stack.

col = '';
if isempty(replicate_table) || ~istable(replicate_table)
    return
end

candidates = {'subject_id','sub','subject','subjectid','subj','sid'};
vars = replicate_table.Properties.VariableNames;
vars_lower = lower(vars);
for i = 1:numel(candidates)
    hit = find(strcmp(vars_lower, candidates{i}), 1, 'first');
    if ~isempty(hit)
        col = vars{hit};
        return
    end
end

end
