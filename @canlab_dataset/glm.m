function out = glm(D, Yvarname, Xvarnames, wh_keep)
%
% predict Y from X using GLM
%
% Usage:
% ----------------------------------------------------------------------------------
% out = glm(D, Yvarname, Xvarnames, wh_keep)
% 
% Examples:
% ----------------------------------------------------------------------------------
% out = glm(D, 'DeltaDon_avg', prednames, wh_keep)
%
% Inputs:
% ----------------------------------------------------------------------------------
% D             a canlab_dataset object
% Yvarname      the name of a variable to predict. must be subject level
% Xvarnames     the name(s) of predictor variables. if multiple, put in
%               cell array. must be subject_level
% wh_keep       a vector of 1/0 values to use as wh_keep
%
% Outputs:
% ----------------------------------------------------------------------------------
% same as for glmfit()
%
% % Copyright Tor Wager, 2013

if nargin < 4 || isempty(wh_keep)
    wh_keep = true(size(D.Subj_Level.id)); %everyone
end


[y, ally, levelY] = get_var(D, Yvarname, wh_keep);
[X, allX, levelX] = get_var(D, Xvarnames, wh_keep);

% Print var name(s)
fprintf('Y (outcome): %s\n', Yvarname);
fprintf('X (predictors): ')
fprintf('%s\t', Xvarnames{:});
fprintf('\n')

if levelY == 2 && levelX == 2
    % Multi-level
    out = igls_multicond(ally, allX, 'iter', 10);
    
elseif levelY ~= 1 || levelX ~= 1, error('Vars must be subject level');
    
else
    [b, dev, stat] = glmfit(X, y);
    glm_table(stat, Xvarnames);
    
    out.b = b;
    out.dev = dev;
    out.stat = stat;
    
end


end
