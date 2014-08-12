function holdout_sets = xval_select_holdout_set_categoricalcovs(covs)
    
u = unique(covs, 'rows');
nu = size(u, 1);  % number of unique cells; also size of holdout sets

[N, k] = size(covs);

onevec = ones(N, 1);

% indices for each cell
for i = 1:nu
    cell_indices{i} = all(covs - u(i * onevec, :) == 0, 2);
end

nobs_per_cell = sum(cat(2, cell_indices{:}));
nfolds = max(nobs_per_cell);  % number of folds; max() leaves some unbalanced holdout sets; min keeps balance but does not test all obs

% random selection of one obs per cell for each fold
for i = 1:nu
    wh_obs{i} = find(cell_indices{i}); 
    wh_obs{i} = wh_obs{i}(randperm(nobs_per_cell(i)));
end

% deal observations into folds
for i = 1:nfolds
    wh = [];
    for j = 1:nu
        if nobs_per_cell(j) >= i 
            wh = [wh wh_obs{j}(i)];
        end
    end
    
    holdout_sets{i} = logical(false * onevec);
    holdout_sets{i}(wh) = true;
    
end

end  % function

% checking stuff
% allsets = cat(2, holdout_set{:});
% sum(allsets)
% sum(allsets, 2)
% covs(holdout_set{1}, :)
% covs(holdout_set{2}, :)
% covs(holdout_set{3}, :)
% covs(holdout_set{4}, :)
% covs(holdout_set{end}, :)