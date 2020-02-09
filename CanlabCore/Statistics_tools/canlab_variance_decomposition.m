function OUT = canlab_variance_decomposition(y, X, X_type, varargin)
% Decompose and plot model variance explained, including categorical (e.g., participant) and continuous variables
%
% OUT = canlab_variance_decomposition(y, X, X_type, varargin)
%
% Example:
% load carsmall
% y = MPG;
% X = [Model_Year Weight];
% X_type = {'categorical' 'continuous'};
% OUT = canlab_variance_decomposition(y, X, X_type);
%
% colors = cell vector with {r g b} for vars 1...k, then shared, then error

% Inputs and names
% -------------------------------------------------------------------------

doverbose = true;
doplots = true;
varnames = {};
for i = 1:size(X, 2), varnames{i} = sprintf('Var_%d', i); end

colors = scn_standard_colors(size(X, 2) + 1);
colors{end + 1} = [.5 .5 .5]; % for error

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

[resid_var_full_model, ~, dfe, df_model] = get_residual_variance(yn, full_xbasis);

var_explained_full_model = 1 - resid_var_full_model / var_y;



% Loop through variables and get resid var, total variance explained for each
for i = 1:size(X, 2)
    
    % Total variance explained by this variable (i) alone
    % -------------------------------------------------------------------------
    
    % Get residual variance after removing this variable (i) only
    resid_var(i) = get_residual_variance(yn, xbasis{i});
    
    total_var_explained(i) = 1 - resid_var(i) / var_y;
    
    % Unique variance explained
    % -------------------------------------------------------------------------
    
    % Get partial residuals after removing ALL OTHER variables
    
    xbasis_other_vars = xbasis;
    xbasis_other_vars(i) = [];                          % remove target variable
    xbasis_other_vars = cat(2, xbasis_other_vars{:});   % Concatenate into design matrix
    
    resid_var_without_var_i = get_residual_variance(yn, xbasis_other_vars);
    
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
    
    create_figure('pie charts')
    
    subplot(1, 2, 1);
    
    % Pie chart of total variance excluding participant-level intercepts
    piedata = double([unique_var_explained shared_var_explained 1-var_explained_full_model]);
    h = wani_pie(piedata, 'hole', 'colors', colors, varargin{:});
    
    title('% of total variance');
    
    subplot(1, 2, 2);
    
    % Pie chart of all explained variance excluding participant-level intercepts
    piedata = double([unique_var_explained shared_var_explained] ./ var_explained_full_model);
    h = wani_pie(piedata, 'hole', 'colors', colors, varargin{:});
    
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



function [resid_var, x_resid, dfe, k] = get_residual_variance(y, xbasis)
% y = outcome
% xbasis = set of one or more regressors for this variable
% assumes all NaNs removed

xbasis = double(xbasis);

xbasis(:, end+1) = 1;                       % Add intercept

x_resid = y - xbasis * pinv(xbasis) * y;  % Residual variance after removing regressor(s)

[n, k] = size(xbasis);                      % Non-Nan obs (n) and k regressor(s) for this variable/factor
% k = 2 (intercept + regressor) unless this is a
% categorical factor

dfe = n - k;                                % intercept accounted for here

resid_var = get_variance(x_resid, dfe);   % Residual uUnexplained) variance after removing xbasis

end



% Define function to get variance; like var() but can input dfe depending on model used to remove other effects
function var_x = get_variance(x, dfe)

msefun = @(x, dfe) x' * x / dfe;
x = double(x - mean(x));                       % Remove mean from each column

var_x = msefun(x, dfe);

end

