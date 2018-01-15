function [X, e] = create_block_design(scanLength, TR, nconditions, blockduration, HPlength, dononlin)
% [X, e] = create_block_design(scanLength, TR, nconditions, blockduration, HPlength, dononlin)
%
% alternating blocks of num_conditions nconditions
%
% scanLength in sec
% enter num_conditions in conditions (e.g. 2)
% enter block_len in sec
%
% X = high-pass filtered design matrix
% e = efficiency
%
% see plotDesign onsets2fmridesign  create_orthogonal_contrast_set


% ----------------------------------------------------------------
% * Fixed parameters
% ----------------------------------------------------------------
%scanLength = 400; % in sec
dohpfilt = 1;

%len = ceil(scanLength .* TR);
ons = cell(1, nconditions);

if dononlin
    nonlinstr = 'nonlinsaturation';
else
    nonlinstr = 'nononlin';
end

% Block onsets
% ----------------------------------------------------------------

for i = 1:nconditions
    
    startval = (i - 1) .* blockduration;
    ons{i} = [startval:(nconditions .* blockduration + 1):scanLength]';  % + 1: 1 sec to start next
    
    ons{i}(:, 2) = blockduration;
    
end

% contrasts: differences among block types
% ----------------------------------------------------------------

contrasts = create_orthogonal_contrast_set(nconditions);
contrasts(:, end+1) = 0; % for intercept

% there should be one column per condition in your design, *including* the
% intercept. One row per contrast.

    X = onsets2fmridesign(ons, TR, scanLength, 'hrf', nonlinstr);


    % high-pass filtering
    if dohpfilt
        % filter all columns except intercept
        scanLength_in_TRs = size(X, 1);
        X(:, 1:end-1) = hpfilter(X(:, 1:end-1), TR, HPlength, scanLength_in_TRs);
    end
    
    plotDesign(ons,[], TR, 'samefig', 'durs', 1, nonlinstr);
    set(gca, 'XLim', [0 scanLength], 'XTick', round(linspace(0, scanLength, 10)));
    
    % test efficiency
    
    e = calcEfficiency(ones(1, size(contrasts, 1)), contrasts, pinv(X), []);
    
    % related to (same as without contrasts, autocorr, filtering):
    % 1 ./ diag(inv(X' * X))

    
end % function




