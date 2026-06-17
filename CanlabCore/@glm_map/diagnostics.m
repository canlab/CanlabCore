function obj = diagnostics(obj, varargin)
% Compute and report design diagnostics for a glm_map object.
%
% Evaluates the conditioning of the design matrix X (and contrasts C):
% variance inflation factors (VIF) per regressor, contrast VIFs (cVIF),
% per-observation leverage, condition number, rank deficiency, and a
% redundant/near-collinear row report. Results are stored back into the
% object and (optionally) printed.
%
% :Usage:
% ::
%
%     obj = diagnostics(obj, varargin)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design matrix available (obj.X non-empty).
%
% :Optional Inputs:
%
%   **'doverbose' / 'noverbose':**
%        Print a summary table (default true).
%
%   **'vif_thresh', [value]:**
%        VIF warning threshold (default 4).
%
% :Outputs:
%
%   **obj:**
%        glm_map with vif, contrast_vif, leverages, condition_number,
%        rank_deficient, collinearity_report, and warnings populated.
%
% :Examples:
% ::
%
%     g = glm_map('X', [ones(30,1) zscore((1:30)')], 'level', 2);
%     g = diagnostics(g);
%
% :See also:
%   - VIF, cVIF, fmri_data.regress
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation. Operates on the design only (no fit
%    required), so it can be used to screen a design before fitting.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
doverbose = true;
if any(strcmpi(varargin, 'noverbose')), doverbose = false; end

vif_thresh = 4;
wh = find(strcmpi(varargin, 'vif_thresh'));
if ~isempty(wh), vif_thresh = varargin{wh(1) + 1}; end

X = obj.X;
if isempty(X)
    error('glm_map:NoDesign', 'No design matrix available (obj.X is empty). Build or supply a design first.');
end

mywarnings = {};

% -------------------------------------------------------------------------
% Variance inflation factors (per regressor)
% -------------------------------------------------------------------------
obj.vif = VIF(X);

if any(obj.vif > vif_thresh)
    mywarnings{end + 1} = sprintf(['Design multicollinearity: %d regressor(s) have VIF > %g. ' ...
        'Check obj.vif and obj.regressor_names.'], sum(obj.vif > vif_thresh), vif_thresh);
end

% -------------------------------------------------------------------------
% Contrast variance inflation factors (per contrast), if contrasts defined
% -------------------------------------------------------------------------
if ~isempty(obj.contrasts)
    if size(obj.contrasts, 1) ~= size(X, 2)
        mywarnings{end + 1} = sprintf(['Contrast matrix has %d rows but design has %d regressors; ' ...
            'skipping contrast VIFs.'], size(obj.contrasts, 1), size(X, 2));
    else
        % cVIF expects one contrast per row -> transpose [P x K] to [K x P]
        obj.contrast_vif = cVIF(X, obj.contrasts');
    end
end

% -------------------------------------------------------------------------
% Leverage (per observation)
% -------------------------------------------------------------------------
H = X * pinv(X);
obj.leverages = diag(H)';

if any(abs(zscore(obj.leverages)) >= 3)
    mywarnings{end + 1} = ['Some observations have high leverage (abs(z(leverage)) >= 3); ' ...
        'the fit may be unstable. Check obj.leverages.'];
end

% -------------------------------------------------------------------------
% Conditioning / rank
% -------------------------------------------------------------------------
obj.condition_number = cond(X);
obj.rank_deficient   = rank(X) < size(X, 2);

if obj.rank_deficient
    mywarnings{end + 1} = 'Design matrix X is rank deficient (rank(X) < number of regressors).';
end

% -------------------------------------------------------------------------
% Redundant / near-collinear column report
% -------------------------------------------------------------------------
report = struct();
report.vif_threshold      = vif_thresh;
report.high_vif_columns   = find(obj.vif > vif_thresh);

% Duplicate (identical) columns
ncol = size(X, 2);
dup_pairs = [];
for a = 1:ncol - 1
    for b = a + 1:ncol
        if isequal(X(:, a), X(:, b))
            dup_pairs(end + 1, :) = [a b]; %#ok<AGROW>
        end
    end
end
report.duplicate_column_pairs = dup_pairs;
if ~isempty(dup_pairs)
    mywarnings{end + 1} = sprintf('%d pair(s) of identical design columns detected (see obj.collinearity_report).', size(dup_pairs, 1));
end

% Near-collinear pairs by pairwise correlation magnitude
R = corrcoef(X);
R(logical(eye(ncol))) = 0;
[ia, ib] = find(triu(abs(R) > 0.95, 1));
report.high_correlation_pairs = [ia ib];

obj.collinearity_report = report;

% -------------------------------------------------------------------------
% Store warnings
% -------------------------------------------------------------------------
obj.warnings = [obj.warnings(:); mywarnings(:)]';
obj.history{end + 1} = 'diagnostics: computed VIF, cVIF, leverage, condition number, collinearity report';

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
if doverbose
    rn = obj.regressor_names;
    fprintf('\n  glm_map diagnostics\n  %s\n', repmat('-', 1, 50));
    fprintf('  %-28s %s\n', 'Regressor', 'VIF');
    for i = 1:numel(obj.vif)
        if i <= numel(rn) && ~isempty(rn{i}), name = rn{i}; else, name = sprintf('R%d', i); end
        flag = ''; if obj.vif(i) > vif_thresh, flag = '  <-- high'; end
        fprintf('  %-28s %6.2f%s\n', name, obj.vif(i), flag);
    end
    if ~isempty(obj.contrast_vif)
        fprintf('  %s\n  Contrast VIFs:\n', repmat('-', 1, 50));
        for i = 1:numel(obj.contrast_vif)
            if i <= numel(obj.contrast_names) && ~isempty(obj.contrast_names{i})
                name = obj.contrast_names{i};
            else
                name = sprintf('Con%d', i);
            end
            fprintf('  %-28s %6.2f\n', name, obj.contrast_vif(i));
        end
    end
    fprintf('  %s\n', repmat('-', 1, 50));
    fprintf('  condition number : %.2f\n', obj.condition_number);
    fprintf('  rank deficient   : %d\n', obj.rank_deficient);
    if ~isempty(mywarnings)
        fprintf('  %d warning(s):\n', numel(mywarnings));
        for i = 1:numel(mywarnings), fprintf('    - %s\n', mywarnings{i}); end
    end
    fprintf('\n');
end

end % diagnostics
