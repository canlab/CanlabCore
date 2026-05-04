function Runc = hrf_build_stick_functions(E, cond_names, TR, n_tp)
%HRF_BUILD_STICK_FUNCTIONS Build one stick function per condition.
Runc = cell(1, numel(cond_names));
for c = 1:numel(cond_names)
    run = zeros(n_tp, 1);
    idx = strcmp(E.trial_type, cond_names{c});
    onsets = E.onset(idx);
    onset_idx = round(onsets ./ TR) + 1;
    onset_idx = onset_idx(onset_idx >= 1 & onset_idx <= n_tp);
    run(onset_idx) = 1;
    Runc{c} = run;
end
end
