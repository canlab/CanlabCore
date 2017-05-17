function [Xscaled, Xcell_scaled] = scale_within_condition(X, conditioning_var, varargin)
% Scales a variable by independently scaling subsets of observations, e.g.,
% within a person, study, or condition. The input conditioning_var
% specifies the subsets of observations on which to operate.
%
% :Usage:
% ::
%
%     [Xscaled, Xcell_scaled] = scale_within_condition(X, conditioning_var, [optional inputs]) 
%
%     - omits NaNs and replaces them after scaling
%     - you can input any anonymous function handle for scaling/operation
%     - (can do any operation on subsets of elements defined by conditioning_var)
% ..
%     Author and copyright information:
%
%     Copyright (C) 2017 Tor Wager
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
%
% :Inputs:
%
%   **X:**
%        a variable, n x 1, or possibly a matrix if your scaling function
%        operates correctly on matrices.
%
% **conditioning_var:**
%       A variable of size n x 1 that contains integer codes (e.g., 1...2...3...) 
%       for subsets of observations belonging to the same group (e.g., subject, study, condition) 
%       The scaling operation will be performed within these subsets.
%
% :Optional Inputs:
%
%   **scaling function:**
%        anonymous function specifying your scaling function, to be
%        performed on the input matrix X.
%        e.g., @zscore
%
%   **zscore**
%       Keyword to use z-score scaling
%
% Examples:
% 
% X = normrnd(0, 1, 99, 1);
% C = [ones(33,1); 2*ones(33, 1); 3*ones(33, 1)];
% Xc = scale_within_condition(X, C); % Norm rank: default
% Xz = scale_within_condition(X, C, 'zscore'); % Z-score
% Xz2 = scale_within_condition(X, C, @zscore); % Z-score
% Xcentered = scale_within_condition(X, C, @(X) scale(X, 1)); % Center X within levels of C
% Xmad = scale_within_condition(X, C, @(X) X ./ mad(X)); % Divide by median abs deviation within levels of C


% ..
%    DEFAULTS AND INPUTS
% ..

% "Nonparametric" : Normalized rank data within-study
% Divide by length (N observations) to keep range same for all studies
scalefun = @(X) rankdata(X) ./ length(X);  % default: Normalized rank

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            %case 'zscore', rowsz = varargin{i+1}; varargin{i+1} = [];
            case 'zscore', scalefun = @(X) zscore(X);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
        
    elseif isa(varargin{i}, 'function_handle')
        
        scalefun = varargin{i};
        
    end
end

% get unique values of conditioning variable
if ~iscolumn(conditioning_var)
    conditioning_var = conditioning_var';
end

u = unique(conditioning_var)';
n_u = length(u);
n_x = size(X, 2);

Xscaled = NaN .* ones(size(X));

for i = 1:n_u
    
    % Select cases
    wh = conditioning_var == u(i);
    
    Xcell_scaled{i} = X(wh, :);
    
    for j = 1:n_x
        % for each column of X
        % Remove NaNs
        [wasnan{i}{j}, xtmp] = nanremove(Xcell_scaled{i}(:, j));
        
        xtmp = scalefun(xtmp);
        
        Xcell_scaled{i}(:, j) = naninsert(wasnan{i}{j}, xtmp);
        
    end % columns
    
    Xscaled(wh, :) = Xcell_scaled{i};
end



end % function