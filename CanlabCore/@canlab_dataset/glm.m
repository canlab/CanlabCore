function out = glm(D, Yvarname, Xvarnames, wh_keep)
% predict Y from X using GLM
%
% :Usage:
% ::
%
%    out = glm(D, Yvarname, Xvarnames, wh_keep)
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
%        the name of a variable to predict. must be subject level
%
%   **Xvarnames:**
%        the name(s) of predictor variables. if multiple, put in
%        cell array. must be subject_level
%
%   **wh_keep:**
%        a logical vector of 1/0 values
%
%
% :Outputs:
%   **out:**
%       structure containing same output as for glmfit()
%       out.b: a vector of coefficient estimates
%       out.dev: the deviance of the fit
%       out.stat: see glmfit documentation for stat structure fields
%
%
% :Examples:
% ::
%
%     out = glm(D, 'DeltaDon_avg', prednames, wh_keep)
%

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
