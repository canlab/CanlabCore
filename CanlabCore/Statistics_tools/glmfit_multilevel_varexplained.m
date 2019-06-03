function vardecomp = glmfit_multilevel_varexplained(X, Y, b, varargin)
% Variance explained estimates, table, and pie chart for multi-level linear model
%
% :Usage:
% ::
%
%     vardecomp =  glmfit_multilevel_varexplained(XX,YY, b, ['noplots' 'notable'])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) <year>  <name of author>
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
%   **b:**
%        is a column vector of beta weights for fixed effects (i.e., returned by glmfit_multilevel: stats.beta')
%        MUST be [intercept followed by beta weights]
%
%   **X:**
%        is a cell vector containing the design matrix for each subject, each subject in a cell. 
%        *columns must code for the same variable for all subjects*
%
%   **Y:**
%        is the data vector for each subject in each cell
%
% :Optional Inputs:
%
%   **'noplots', 'notable':**
%        suppress output
%
%   **'colors':**
%        followed by n x 3 vector of rgb colors
%        other inputs to wani_pie OK too.
%
% :Output:
%
%   **vardecomp:**
%        structure of variance explained and description
%
% :See also:
% glmfit_multilevel, igls_multicond, glmfit_general, wani_pie
%
% ..
%    Programmers' notes:
%    Created by Tor Wager and Anjali Krishnan, May 2015
% ..

% ..
%    DEFAULTS AND INPUTS
% ..

doplots = 1;
dotable = 1;

% Varargin
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'noplot', 'noplots'}
                doplots = 0;
                
            case 'notable'
                dotable = 0;
                
            case 'colors'
                % do nothing for compatibility...
            otherwise
                fprintf('Warning! Unknown input string option: %s', varargin{i});  
        end
    end
end


%% -------------------------------------------------------------------------
% SET UP ANALYSIS
% -------------------------------------------------------------------------

% Remove random intercepts per participants by mean-centering Y
% This variance is not relevant for decomposition
Yc = cellfun(@(x) scale(x, 1), Y, 'UniformOutput', false);

% Concatenate all the values across participants
Ycat = cat(1, Yc{:});
Xcat = cat(1, X{:});

% Get number of total observations across participants
nobs = size(Xcat, 1);

% Include intercept for overall model
Xcat = [ones(nobs, 1) Xcat];

% Compute the model fit with betas obtained from glmfit_multilevel
modelfit = Xcat * b;

% ***SHOULD PROBABLY MEAN-CENTER PREDICTED VALUES, BECAUSE WE MEAN-CENTERED
% Y***
% E.G., modelfit = cellfun(@(x) scale(x, 1), Y, 'UniformOutput', false);

 
%% -------------------------------------------------------------------------
% Compute variances components 
% -------------------------------------------------------------------------

% Total variance (without intercept)
vardecomp.vary = var(Ycat);               

% Explained variance (fixed effect variables only, random intercepts removed)
%Even though there is intercept in overall model,  it is a constant (no effect on explained variance)
vardecomp.varexpfixed = var(modelfit);    

% Residual variance including measurement error (sigma) and variance due to random slopes
vardecomp.varresid = vardecomp.vary - vardecomp.varexpfixed;  

% Partial variance: variance associate with each independent variable (I.V.)
% Computed as the product of the I.V. and the corresponding beta across all
% observations and participants
partialfit = Xcat * diag(b);

% Avoid machine precision error in some rare cases (?) by normalizing:
vpf = var(partialfit); 
vpf = vpf .* var(modelfit) ./ sum(vpf);

% Ignore intercept as it has no variance
partialfit = partialfit(:, 2:end); 

% Unique partial variance attributed to each I.V.
vardecomp.varbeta = vpf(2:end);      

% Shared variance between all the I.V.s
% Computed as : varexpfixed = sum(varbeta) + varshared
vardecomp.varshared = vardecomp.varexpfixed - sum(vardecomp.varbeta);

% Check variances, total variance must be sum of all variance components
% vary = varexpfixed + varresid = sum(varbeta) + varshared + varresid

%% -------------------------------------------------------------------------
% PERCENTAGES
% -------------------------------------------------------------------------

% Percentage of total variance (excluding participant-level intercept) attributable to each I.V., and shared variance

vardecomp.tot_varbeta_percent = 100 .* ([vardecomp.varbeta vardecomp.varshared vardecomp.varresid] ./ vardecomp.vary);
vardecomp.tot_varbeta_names = {'1 to n-2 columns = variance for each I.V.' '2nd to last column = shared variance' 'last column = residual variance'};

% Percentage of explained variance (excluding participant-level intercept) attributable to each component

vardecomp.exp_varbeta_percent = 100 .* ([vardecomp.varbeta vardecomp.varshared] ./ vardecomp.varexpfixed);
vardecomp.exp_varbeta_names = {'1 to n-1 columns = variance for each I.V.' 'last column = shared variance'};

%% -------------------------------------------------------------------------
% PLOTS AND OUTPUT
% -------------------------------------------------------------------------

if doplots
    
    create_figure('pie charts')
    
subplot(1, 2, 1);

% Pie chart of total variance excluding participant-level intercepts
h = wani_pie(double(vardecomp.tot_varbeta_percent), 'hole', varargin{:});
title('% of total variance');

subplot(1, 2, 2);

% Pie chart of all explained variance excluding participant-level intercepts
h = wani_pie(double(vardecomp.exp_varbeta_percent), 'hole', varargin{:});

title('% of explained variance');

end

if dotable
    
    fprintf('Variance estimates:\n')
    fprintf('var(y) = var(model) + var(resid)\n')
    fprintf('%3.3f = %3.3f + %3.3f\n\n',vardecomp.vary, vardecomp.varexpfixed, vardecomp.varresid );
    
    x = [vardecomp.varbeta vardecomp.varshared; vardecomp.tot_varbeta_percent(1:end-1); vardecomp.exp_varbeta_percent]';
    
    for i = 1:length(vardecomp.varbeta)
        cnames{i} = sprintf('b%d', i);
    end
    cnames{end+1} = 'shared';
    
    print_matrix(x, {'Var' '% of Total' '% of Model'}, cnames);
    
end


%% In-line Functions
% pie_plot: Makes pie-charts for the variance decomposition and heavily
% borrows from wani_pie written by Wani Woo in April 2015
% NOTE: wani_pie.m is stand-alone function
% 
%     function h = wani_pie(X, varargin)
%         % Draw a little better pie chart
%         %
%         % Usage:
%         % -------------------------------------------------------------------------
%         % h = wani_pie(X, varargin)
%         %
%         % Inputs:
%         % -------------------------------------------------------------------------
%         % X a vector
%         %
%         % Optional inputs: Enter keyword followed by variable with values
%         % 'cols' colors N x 3 {default: using colors from microsoft office}
%         % 'notext' no text for percentage {default: false}
%         % 'fontsize' font size for percentage {default: 15}
%         % 'hole' add a hole in the middle of the pie chart {default: no hole}
%         % 'hole_size' specify the size of the middle hole {default: 5000}
%         % 'outline'
%         % 'outlinecol'
%         % 'outlinewidth'
%         %
%         % Outputs:
%         % -------------------------------------------------------------------------
%         % h graphic handles
%         %
%         % Examples:
%         % -------------------------------------------------------------------------
%         % % data
%         % X = rand(10,1);
%         % h = wani_pie(X, 'notext', 'hole')
%         %
%         % savename = 'example_pie.pdf';
%         %
%         % try
%         % pagesetup(gcf);
%         % saveas(gcf, savename);
%         % catch
%         % pagesetup(gcf);
%         % saveas(gcf, savename);
%         % end
%         %
%         % -------------------------------------------------------------------------
%         % Copyright (C) 2015 Wani Woo
%         cols = [0.0902 0.2157 0.3686
%                 0.2157 0.3765 0.5725
%                 0.5843 0.2157 0.2078
%                 0.4667 0.5765 0.2353
%                 0.3765 0.2902 0.4824
%                 0.1922 0.5216 0.6118
%                 0.8941 0.4235 0.0392];
%             
%         eps = 0.0001;
%         
%         if numel(X) > 7
%             cols = [cols(randperm(7),:);
%                 cols(randperm(7),:);
%                 cols(randperm(7),:)];
%         else
%             cols = cols(randperm(7),:);
%         end
%         
%         dotext = 1;
%         fs = 15;
%         hs = 5000;
%         dohole = 0;
%         doout = 0;
%         outlinecol = [0 0 0];
%         outlinewidth = 1.2;
%         for i = 1:length(varargin)
%             if ischar(varargin{i})
%                 switch varargin{i}
%                     % functional commands
%                     case {'cols'}
%                         cols = varargin{i+1};
%                     case {'notext'}
%                         dotext = 0;
%                     case {'fontsize'}
%                         fs = varargin{i+1};
%                     case {'hole'}
%                         dohole = 1;
%                     case {'hole_size'}
%                         hs = varargin{i+1};
%                     case {'outline'}
%                         doout = 1;
%                     case {'outlinecol'}
%                         outlinecol = varargin{i+1};
%                     case {'outlinewidth'}
%                         outlinewidth = varargin{i+1};
%                 end
%             end
%         end
%         
%         % avoid exactly-zero values, which will give error
%         X(X == 0) = min(X(X ~= 0) ./ 1000); 
% 
%         h = pie(X);
%         set(gcf, 'color', 'w', 'position', [360 393 389 305]);
%         
%         for i = 1:numel(X)
%             set(h(2*i-1), 'facecolor', cols(i,:), 'edgecolor', 'none');
%             hold on;
%             if dotext
%                 set(h(2*i), 'fontSize', fs);
%             else
%                 set(h(2*i), 'String', '');
%             end
%         end
%         
%         if dohole
%             scatter(0, 0, hs, 'w', 'filled');
%         end
%         
%         if doout
%             scatter(0, 0, 43500, outlinecol, 'linewidth', outlinewidth);
%         end
%     end
end
