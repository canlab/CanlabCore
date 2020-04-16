function OUT = canlab_variance_decomposition(y, X, X_type, varargin)
% Decompose and plot model variance explained, including categorical (e.g., participant) and continuous variables
%
% OUT = canlab_variance_decomposition(y, X, X_type, varargin)
%
% :Inputs:
%
%   **param1:**
%        description of param1
%
%   **param2:**
%        description of param2
%
% :Optional Inputs:
%   **'colors':**
%       Followed by cell vector with {r g b} for vars 1...k, then shared, then error
%
%   **'names', 'VarNames':**
%       Followed by variable names. Omit to show
%       percent variance explained instead of variable name
%
%   **'prediction_r2':**
%       Do not re-estimate intercept and betas, and instead calculate
%       absolute model fit between y and X. Use if X is y-hat from a
%       pre-trained model and you want to evaluate absolute prediction r^2
%       Betas will be fixed at 1, and the model has 0 model df.
%
%
% :Outputs:
%
%   **OUT:**
%        A structure with various variance estimates, and a table with
%        total and unique variance explained by each column (predictor) in X.
%
% Notes: Using canlab_variance_decomposition to assess cross-validated or
% pre-trained model prediction
% -------------------------------------------------------------------------
% Two options:
% 1. Pass in y = observed, X = cross-validated or test-sample model
% predictions (y-hat from one or more models). Multiple columns of X would
% be used to test multiple, potentially correlated model predictions from
% different models simultaneously.
%
% By default, canlab_variance_decomposition re-fits a model (X) to y, so it?s 
% not exactly the prediction R^2 in Schienost 2019 (Rule 5). That is, it will 
% fit y ~ X and use the fits to get variance explained for each. If X is a 
% cross-validated or test-sample model prediction, this function will re-estimate 
% the slope and intercept, which will allow rescaling between the observed(y) and 
% predicted (y-hat/X). This may be useful if the scale is different between 
% training and test and you want to allow this minor flexibility (at the expense of 
% having positively biased estimates of R^2 which is *very bad*). But if you want 
% to assess prediction R^2 with the exact same parameters as in the original 
% training model, we use the standard approach of using the model fit on training data. 
% In this case, use 'prediction_r2'. Variance explained can be negative.
%
% 2. Pass in y = y-hat from pre-trained model, and X variables that attempt
% to explain the model predictions (y-hat) as a function of other
% variables (e.g., explanatory variables or confounds). Here re-fitting a model is 
% generally approprate, as the scale of y-hat and predictors may not match.
%
% Examples:
% -------------------------------------------------------------------------
%
% load carsmall
% y = MPG;
% X = [Model_Year Weight];
% X_type = {'categorical' 'continuous'};
% OUT = canlab_variance_decomposition(y, X, X_type);
%
     
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Tor Wager
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
% Contributors: Tor Wager, Yoni Ashar, Phil Kragel

% Inputs and names
% -------------------------------------------------------------------------

doverbose = true;
doplots = true;
varnames = {};

colors = scn_standard_colors(size(X, 2) + 1);
colors{end + 1} = [.5 .5 .5]; % for error

% parse input
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case {'names', 'varnames'}
                varnames = varargin{i + 1};
                
            case {'colors', 'color'}
                colors = varargin{i + 1};
                
            case 'prediction_r2'
                % do nothing here; this is passed in varargin{:} to get_residual_variance
                % if it is, betas will be fixed at 1, and the model has 0 model df.
                
            case 'noplots'
                
                doplots=false;
                
            case 'noverbose'
                doverbose = false;
        end
    end
end

% if varnames was not passed in
if isempty(varnames)
    for i = 1:size(X, 2), varnames{i} = sprintf('Var_%d', i); end
end

% Check colors
validateattributes(colors,{'cell'},{@(x) length(x) >= size(X, 2) + 2});
    
% X_type must be ROW vector
if iscolumn(X_type), X_type = X_type'; end



% Remove NaNs casewise from ALL regressors, for full model

[wasnan, Xn, yn] = nanremove(X, y);    % Remove NaNs casewise from joint X, y

% Get total variance in y
% -------------------------------------------------------------------------

var_y = get_variance(yn, size(yn, 1) - 1); % variation around mean; MST; same as nanvar(y)

% Get regressors
% -------------------------------------------------------------------------

% For each variable entered, check if categorical, and expand into an
% orthogonal basis set if so. Keep track of original variables in cells, so
% we can control for other variables (even if categorical) when calculating
% unique variance explained.

xbasis = expand_X_variables(Xn, X_type);

% Full model variance explained
% -------------------------------------------------------------------------
% Needed for calculating unique var explained for each variable

full_xbasis = cat(2, xbasis{:});   % Concatenate into design matrix

[resid_var_full_model, ~, dfe, df_model] = get_residual_variance(yn, full_xbasis, varargin{:});

var_explained_full_model = 1 - resid_var_full_model / var_y;



% Loop through variables and get resid var, total variance explained for each
for i = 1:size(X, 2)
    
    % Total variance explained by this variable (i) alone
    % -------------------------------------------------------------------------
    
    % Get residual variance after removing this variable (i) only
    resid_var(i) = get_residual_variance(yn, xbasis{i}, varargin{:});
    
    total_var_explained(i) = 1 - resid_var(i) / var_y;
    
    % Unique variance explained
    % -------------------------------------------------------------------------
    
    % Get partial residuals after removing ALL OTHER variables
    
    xbasis_other_vars = xbasis;
    xbasis_other_vars(i) = [];                          % remove target variable
    xbasis_other_vars = cat(2, xbasis_other_vars{:});   % Concatenate into design matrix
    
    if  any(strcmp(varargin, 'prediction_r2')) && isempty(xbasis_other_vars)
    resid_var_without_var_i =   var_y; %in this case there is only one predictor...  
    else
    resid_var_without_var_i = get_residual_variance(yn, xbasis_other_vars, varargin{:});
    end
    % Unique variance explained by target variable (i) after removing all other variables
    
    unique_var_explained(i) = (resid_var_without_var_i - resid_var_full_model) ./ var_y;
    
    
end % Loop through variables

% Variance shared across regressors, not unique to any
% -------------------------------------------------------------------------
shared_var_explained = var_explained_full_model - sum(unique_var_explained);

% Define and populate output structure
% -------------------------------------------------------------------------

OUT = struct();

OUT.summary{1} = sprintf('Total y variance = %3.4f, full-model residual variance = %3.4f', var_y, resid_var_full_model);
OUT.summary{2} = sprintf('Full model variance explained: %3.4f', var_explained_full_model);
OUT.summary{3} = sprintf('Variance shared across variables and not explained uniquely: %3.4f', shared_var_explained);

OUT.inputs.X = X;
OUT.inputs.X_type = X_type;
OUT.inputs.y = y;
OUT.inputs.wasnan = wasnan;

OUT.full_model.xbasis = xbasis;
OUT.full_model.xbasis_descrip = 'Cell array of basis regressors for each input variable';

OUT.full_model.var_y = var_y;
OUT.full_model.resid_var_full_model = resid_var_full_model;
OUT.full_model.var_explained_full_model = var_explained_full_model;
OUT.full_model.dfe = dfe;
OUT.full_model.df_model = df_model;

vtable = table(X_type', total_var_explained', unique_var_explained', 'VariableNames', {'Type' 'Total' 'Unique'});
vtable.Properties.RowNames = varnames;

OUT.variables.var_explained_table = vtable;

OUT.variables.total_var_explained = total_var_explained;
OUT.variables.unique_var_explained = unique_var_explained;
OUT.variables.shared_var_explained = shared_var_explained;

% -------------------------------------------------------------------------
% DISPLAY PLOTS AND OUTPUT
% -------------------------------------------------------------------------

if doverbose
    % Print report
    
    canlab_print_legend_text(OUT.summary{:});
    
    disp('Variance explained by each variable:');
    disp(OUT.variables.var_explained_table);
    
end

if doplots
    
    create_figure('pie charts');
    
    subplot(1, 2, 1);
    
    plotlabels = format_strings_for_legend({varnames{:} 'Shared' 'Unexplained'});
    
    % Pie chart of total variance excluding participant-level intercepts
    piedata = double([unique_var_explained shared_var_explained 1-var_explained_full_model]);
    h = wani_pie(piedata, 'hole', 'colors', colors, 'labels', plotlabels, varargin{:});
    
    title('% of total variance');
    
    subplot(1, 2, 2);
    
    % Pie chart of all explained variance excluding participant-level intercepts
    piedata = double([unique_var_explained shared_var_explained] ./ var_explained_full_model);
    %h = wani_pie(piedata, 'hole', 'colors', colors, 'labels', {varnames{:} 'Shared'}, varargin{:});
    h = wani_pie(piedata, 'hole', 'colors', colors); %, 'labels', {varnames{:} 'Shared'}, varargin{:});
    
    title('% of explained variance');
    
end


end % main function


function xbasis = expand_X_variables(X, X_type)
% For each variable entered, check if categorical, and expand into an
% orthogonal basis set if so. Keep track of original variables in cells, so
% we can control for other variables (even if categorical) when calculating
% unique variance explained.
% Assumes NaNs already removed; pass in Xn in main function

xbasis = cell(1, size(X, 2));

for i = 1:size(X, 2)
    
    x = X(:, i);
    x_type = X_type{i}; % 'continuous' or 'categorical'
    
    
    % Create orthgonal basis set spanning space of levels of categorical factor
    
    switch lower(x_type)
        case 'categorical'
            xbasis{i} = condf2basisset(x);
            
        case 'continuous'
            % do nothing - already a continuous regressor
            xbasis{i} = x;
            
        otherwise
            error('x_type must be categorical or continuous');
    end
    
end % loop

end % subfunction



function x = condf2basisset(condf)
% Create orthgonal basis set spanning space of levels of categorical factor
%
% x = condf2basisset(condf)
%
% condf is a condition function, with integer codes indicating levels of a factor
% condf is assumed to be a nominal variable; i.e., levels are not ordered
% e.g., [1 1 1 2 2 2 3 3 3]' for 9 observations with 3 levels
%
% x is an orthonormal basis set, an orthogonal set of regressors spanning the space of all
% differences among levels of the factor. Regressors are scaled so that
% their norm (vector length) is 1.
% This is useful for adding factors with nominal levels into a regression
% model, i.e., ANOVA as regression.

indic = condf2indic(condf);
contrasts = create_orthogonal_contrast_set(size(indic, 2));
x = orth(double(indic * contrasts'));

end



function [resid_var, x_resid, dfe, k] = get_residual_variance(y, xbasis, varargin)
% y = outcome
% xbasis = set of one or more regressors for this variable
% assumes all NaNs removed

xbasis = double(xbasis);

if ~isempty(varargin) && any(strcmp(varargin, 'prediction_r2'))
    % Fix all betas to 1
    
    x_resid = y - xbasis * ones(size(xbasis, 2), 1);  % Residual variance after removing regressor(s)

    dfe = size(xbasis, 1);
    k = 0;
    
    
else
    % Fit model to estimate intercept and betas
    
    xbasis(:, end+1) = 1;                       % Add intercept

    x_resid = y - xbasis * pinv(xbasis) * y;  % Residual variance after removing regressor(s)

    [n, k] = size(xbasis);                      % Non-Nan obs (n) and k regressor(s) for this variable/factor
    % k = 2 (intercept + regressor) unless this is a
    % categorical factor

    dfe = n - k;                                % intercept accounted for here

end

resid_var = get_variance(x_resid, dfe);   % Residual uUnexplained) variance after removing xbasis

end



% Define function to get variance; like var() but can input dfe depending on model used to remove other effects
function var_x = get_variance(x, dfe)

msefun = @(x, dfe) x' * x / dfe;
x = double(x - mean(x));                       % Remove mean from each column

var_x = msefun(x, dfe);

end

