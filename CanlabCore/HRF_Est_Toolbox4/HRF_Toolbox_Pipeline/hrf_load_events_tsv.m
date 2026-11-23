function E = hrf_load_events_tsv(events_tsv)
%HRF_LOAD_EVENTS_TSV Read BIDS events.tsv and validate required columns.
T = readtable(events_tsv, 'FileType', 'text', 'Delimiter', '\t');
req = {'onset', 'duration', 'trial_type'};
for i = 1:numel(req)
    if ~ismember(req{i}, T.Properties.VariableNames)
        error('events.tsv is missing required column: %s', req{i});
    end
end
E = table();
E.onset = T.onset;
E.duration = T.duration;
E.trial_type = cellstr(string(T.trial_type));
end
