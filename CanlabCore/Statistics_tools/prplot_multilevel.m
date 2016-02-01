function [X_resid, Y_resid, handles] = prplot_multilevel(Y, X, wh_col)
% Partial correlation plot for multi-level analysis
%
% :Usage:
% ::
%
%     [X_resid, Y_resid, handles] = prplot_multilevel(Y, X, wh_col)
%
% Uses unweighted estimates.
%
% do not enter intercept in X

N = length(X);

for i = 1:N
    
    X_other{i} = X{i};
    X_other{i}(:, wh_col) = [];
    
    X_part{i} = X{i}(:, wh_col);
    
end

wh = true(size(X{1}, 2) + 1, 1); % assume X is same size; k vars + intercept

wh(wh_col + 1) = false;  % omit from betas; consider intercept

stats = glmfit_multilevel(Y, X_other, [], 'noverbose');

for i = 1:N
    
    Y_resid{i} = Y{i} - [ones(size(X_other{i}, 1), 1) X_other{i}] * stats.first_level.beta(:, i);
    
end

stats = glmfit_multilevel(X_part, X_other, [], 'noverbose');

for i = 1:N
    
    X_resid{i} = X_part{i} - [ones(size(X_other{i}, 1), 1) X_other{i}] * stats.first_level.beta(:, i);
    
end


handles.overall = plot(cat(1, X_resid{:}), cat(1, Y_resid{:}), 'ko');

% partial fit plot
for i = 1:N
    
    handles.indiv = plot(X_resid{i}, Y_resid{i},'Color', rand(1, 3));
    
end

handles.refline = refline;
set(handles.refline, 'LineWidth', 3);

end % function
