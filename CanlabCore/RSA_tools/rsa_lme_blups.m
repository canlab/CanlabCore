function T = rsa_lme_blups(mdl)
% rsa_lme_blups  Extract per-subject BLUPs from a fitted LinearMixedModel.
%
% Returns Best Linear Unbiased Predictors (BLUPs) of the random effects per
% subject. For an LME with `(1 | Sub)`, each subject gets their own
% deviation from the population intercept. For `(X | Sub)`, also their own
% deviation from the population slope on X.
%
% Replicates `08072024 Run-Level RDM Analysis with RSA Toolbox.mlx`
% lines 2226-2280.
%
% Usage
% -----
%   mdl = dat.rsa_lme('Y ~ SameCondition + (SameCondition | Sub)', ...);
%   T = rsa_lme_blups(mdl);
%   % T columns: Group | Subject | Term | Estimate | SEPred (if available)
%
% Inputs
% ------
%   mdl  LinearMixedModel
%
% Output
% ------
%   T   long-format table:
%         Group    -- random-effects grouping variable name (e.g. 'Sub')
%         Subject  -- per-group level (subject ID)
%         Term     -- random-effects term name (e.g. '(Intercept)', 'SameCondition')
%         Estimate -- BLUP value for this (subject, term)
%         SEPred   -- standard error of the prediction (when MATLAB provides it)

assert(isa(mdl, 'LinearMixedModel'), 'rsa_lme_blups: requires a LinearMixedModel input.');

% Use the dataset overload of randomEffects to get a table with metadata
[blups, blup_names, blup_stats] = randomEffects(mdl);

% blup_stats is a dataset/table with Group, Level, Name columns describing
% each entry of blups (one row per (group, level, term))
if istable(blup_stats)
    stats_tbl = blup_stats;
else
    stats_tbl = dataset2table(blup_stats);
end

n = numel(blups);

% Coerce to a uniform set of columns
G = stats_tbl.Group;
L = stats_tbl.Level;
N = stats_tbl.Name;

if iscategorical(G), G = cellstr(G); end
if iscategorical(L), L = cellstr(L); end
if iscategorical(N), N = cellstr(N); end

SE = nan(n, 1);
if ismember('SEPred', stats_tbl.Properties.VariableNames)
    SE = double(stats_tbl.SEPred);
end

T = table(G, L, N, blups, SE, ...
    'VariableNames', {'Group', 'Subject', 'Term', 'Estimate', 'SEPred'});

end
