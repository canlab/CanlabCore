function [b, dev, stat] = glm_multilevel(D, Yvarname, Xvarnames, wh_keep)
%
% predict Y from X using GLM
%
% Usage:
% ----------------------------------------------------------------------------------
% [b, dev, stat] = glm(D, 'DeltaDon_avg', prednames, wh_keep)
%
% Inputs:
% ----------------------------------------------------------------------------------
% D             a canlab_dataset object
% Yvarname      the name of a variable to predict. must be event level
% Xvarnames     the name(s) of predictor variables. if multiple, put in
%               cell array. must be event level
% wh_keep       a vector of 1/0 values to use as wh_keep
%
% Outputs:
% ----------------------------------------------------------------------------------
% 
%
% % Copyright Tor Wager, 2013

[Y, ~, levelY] = get_var(D, Yvarname, wh_keep);
Y = Y';

%% MUST IMPLEMENT GET_VAR FOR MULTIPLE VARS AT AN EVENT LEVEL.  RETURNS A CELL ARRAY FOR EACH PERSON.  DAT BECOMES WARNING STRING, DATCELL IS OF INTEREST.

[X, ~, levelX] = get_var(D, Xvarnames, wh_keep);      

if levelY ~= 2 || levelX ~= 2, error('Vars must be event level'); end



n=size(Y,2);
X1 = cell(1,n);
X2 = mean(Y)'; % matrix of 2nd level preds

if isstruct(varargin{1}) % the "alternative format" described above
    xstruct = varargin{1};
    fields = varargin{2};
    for i = 1:n % each subject
        clear myX
        for j = 1:length(fields) %all the fields
            myX(:, j) = xstruct.(fields{j}){i};
            X1{i} = myX;
        end
    getvif(X1{i})
    end

else
    for i = 1:n % each subject
        clear myX
        for j = 1:length(varargin)
            myX(:, j) = varargin{j}{i};
            X1{i} = myX;
        end
    end
end


stats = glmfit_multilevel(Y, X1, scale(X2, 1), 'weighted', 'noplots', ...
    'names', {'Intrcpt' names{2:end}}, 'beta_names', {'Avg within-ss relationship' ['Effect of Avg. ' names{1}]});

end