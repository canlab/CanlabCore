function T = hrf_summarize_study(results, subject_ids)
%HRF_SUMMARIZE_STUDY Long-format summary table per subject/condition/model.
rows = {};
for s = 1:numel(results)
    r = results{s};
    models = fieldnames(r.fits);
    for m = 1:numel(models)
        mdl = models{m};
        if ~isfield(r.fits.(mdl), 'param'), continue; end
        P = r.fits.(mdl).param;
        for c = 1:numel(r.conditions)
            rows(end+1, :) = {subject_ids{s}, mdl, r.conditions{c}, P(1,c), P(2,c), P(3,c)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'subject','model','condition','amplitude','time_to_peak','width'});
end
