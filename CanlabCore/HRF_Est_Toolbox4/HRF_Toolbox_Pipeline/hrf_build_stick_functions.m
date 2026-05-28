function [Runc, condition_groups] = hrf_build_stick_functions(E, cond_names, TR, n_tp)
%HRF_BUILD_STICK_FUNCTIONS Build one stick function per condition.
available_conditions = unique(cellstr(string(E.trial_type)), 'stable');
condition_groups = hrf_resolve_condition_patterns(available_conditions, cond_names, 'DefaultMode', 'each');

Runc = cell(1, numel(condition_groups));
trial_type = cellstr(string(E.trial_type));
for c = 1:numel(condition_groups)
    run = zeros(n_tp, 1);
    idx = ismember(trial_type, condition_groups(c).matched_conditions);
    onsets = E.onset(idx);
    onset_idx = round(onsets ./ TR) + 1;
    onset_idx = onset_idx(onset_idx >= 1 & onset_idx <= n_tp);
    run(onset_idx) = 1;
    Runc{c} = run;
end
end
