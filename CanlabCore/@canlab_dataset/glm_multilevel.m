function [b, dev, stat] = glm_multilevel(D, Yvarname, Xvarnames, wh_keep)
% Predict Y from X using GLM
%
% :Usage:
% ::
%
%    [b, dev, stat] = glm_multilevel(D, Yvarname, Xvarnames, wh_keep)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2013 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **D:**
%        a canlab_dataset object
%
%   **Yvarname:**
%        the name of a variable to predict. must be event level
%
%   **Xvarnames:**
%        the name(s) of predictor variables. if multiple, put in
%        cell array. must be event level
%
%   **wh_keep:**
%        a logical vector of 1/0 values
%
% :Outputs:
%   **b:**
%       a vector of coefficient estimates (same as for glmfit())
%
%   **dev:**
%       the deviance of the fit (same as for glmfit())
%
%   **stat:**
%       structure containing stats fields (see glmfit() documentation)
%

[Y, ~, levelY] = get_var(D, Yvarname, wh_keep);
Y = Y';

%% MUST IMPLEMENT GET_VAR FOR MULTIPLE VARS AT AN EVENT LEVEL.  RETURNS A CELL ARRAY FOR EACH PERSON.  DAT BECOMES WARNING STRING, DATCELL IS OF INTEREST.

[X, ~, levelX] = get_var(D, Xvarnames, wh_keep);      

if levelY ~= 2 || levelX ~= 2, error('Vars must be event level'); end



n=size(Y,2);
X1 = cell(1,n);
X2 = mean(Y)'; % matrix of 2nd level preds

%{
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
%}
    for i = 1:n % each subject
        clear myX
        for j = 1:length(varargin)
            myX(:, j) = varargin{j}{i};
            X1{i} = myX;
        end
    end
%end


stats = glmfit_multilevel(Y, X1, scale(X2, 1), 'weighted', 'noplots', ...
    'names', {'Intrcpt' names{2:end}}, 'beta_names', {'Avg within-ss relationship' ['Effect of Avg. ' names{1}]});

end
