function [nuisance_ratio, rsquare_design, rsquare_nuis] = scn_component_rsquare(V, nuisanceX, designX)
% Print a table of r-square values (variance explained) for each of V data
% vectors by nuisance (mvmt, physio) and task-related predictors
%
% Designed to work with components
%
% :Examples:
% ::
%
%     % Typical operation
%     scn_component_rsquare(compscore, movement_params(1:157, :), X(1:157, :));
%
%     % No design
%     scn_component_rsquare(compscore, movement_params(1:157, :));
%
%     % Neither design nor nuisance, uses linear drift
%     scn_component_rsquare(compscore, []);
%
% ..
%    Tor Wager, Feb 2008
% ..

[t, m] = size(V);

if nargin < 3, designX = []; end

if isempty(nuisanceX) 
    disp('Looking for nuisance covariates, but none found. Using linear drift.');
    nuisanceX = (1:t)'; 
end
    
disp('Variance in each component explained:')

% Nuisance
rsquare_nuis = get_rsquare(V, nuisanceX)';

if ~isempty(designX)
rsquare_design = get_rsquare(V, designX)';
else
    rsquare_design = .01 * ones(m, 1);
end

nuisance_ratio = rsquare_nuis ./ rsquare_design;

nuis_related = find(nuisance_ratio > 2);
task_related = find(nuisance_ratio < 1);

%rank_badness = sort(nuisance_ratio, 1, 'descend');
disp('All components');
print_matrix([(1:m)' rsquare_design rsquare_nuis nuisance_ratio], {'Comp.' 'R^2 Task' 'R^2 Nuisance' 'Ratio'});
fprintf('\n');

disp('Most task-related');
if isempty(designX)
    disp('NO DESIGN INFORMATION.');
else
print_matrix([task_related rsquare_design(task_related) rsquare_nuis(task_related) nuisance_ratio(task_related)], {'Comp.' 'R^2 Task' 'R^2 Nuisance' 'Ratio'});
end
fprintf('\n');

disp('Most nuisance-related');
print_matrix([nuis_related rsquare_design(nuis_related) rsquare_nuis(nuis_related) nuisance_ratio(nuis_related)], {'Comp.' 'R^2 Task' 'R^2 Nuisance' 'Ratio'});
fprintf('\n');


%fprintf('\tComp %3.0f : %3.0f%%\n', j, rsquare(j)*100);
end

function rsquare = get_rsquare(V, X)

    [t, m] = size(V);
    
    for j = 1:m
    % Center component to avoid counting intercept in r-square?
    % Don't have to, doesn't matter because var operator is 2nd moment
    y = V(:, j);  
    b =  X \ y;

    fits = X * b;

    rsquare(j) = var(fits) / var(y);
    
    end

    
end
