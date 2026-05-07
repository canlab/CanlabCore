function ax = plot_hrf_study_by_subject(study, varargin)
%PLOT_HRF_STUDY_BY_SUBJECT Plot per-subject HRFs separated by subject.
p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.parse(study, varargin{:});
opts = p.Results;

model_name = char(opts.Model);
figure; ax = axes; hold on;
colors = lines(numel(study.results));
for s = 1:numel(study.results)
    r = study.results{s};
    if isempty(r), continue; end
    if isempty(opts.Condition)
        c = 1;
    else
        c = find(strcmp(r.conditions, char(opts.Condition)), 1);
    end
    if isempty(c), continue; end
    y = r.fits.(model_name).hrf(:, c);
    plot(y, 'Color', colors(s,:), 'DisplayName', study.subject_ids{s});
end
legend('Interpreter','none');
title(sprintf('%s by subject', upper(model_name)));
hline(0,'k-');
end
