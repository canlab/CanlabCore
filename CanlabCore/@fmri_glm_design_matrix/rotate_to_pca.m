function [obj, x_eig, c_eig] = rotate_to_pca(obj)
% Rotate design matrix columns within all conditions to principal component projection.


eigval_cutoff = .1; % 100*eps;  % or 1

meth = 'nuis only';

switch meth
    case 'all'
        
        c = obj.xX.cond_assignments;
        
        % add column for all nuisance covs
        c(:, end + 1) = (~any(c'))';
        
        % here, if startval = size(c, 2), then
        % it will do only nuisance covariates.  This is the default behavior.
        % if startval were 1, it would do ALL conditions.
        startval = 1;
        
    case 'nuis only'
        
        % this for nuisance only instead
        c = false(size(obj.xX.X, 1), 1);
        c(obj.xX.iC) = true;
        
        startval = size(c, 2);
        
end


[c_eig, x_eig] = deal(cell(1, size(c, 2)));



for j = startval:size(c, 2)
    % for each condition
    
    % columns for condition j
    x = obj.xX.X(:, c(:, j));
    
    % for determining which to save and which are empty
    x2 = scale(x);
    
    [v, sc, lam] = princomp(x2, 'econ');
    wh_empty = lam < eigval_cutoff;
    
    v = princomp(x, 'econ');
    
    % eliminate empty columns
    v(:, wh_empty) = [];
    x = x * v;
    
    % re-build c matrix
    c_eig{j} = ones(size(x, 2), 1);
    
    x_eig{j} = x;
    
end

x_eig = cat(2, x_eig{:});
c_eig = blkdiag(c_eig{:});

switch meth
    case 'all'
        
        
        
        
    case 'nuis only'
        wh_to_replace = obj.xX.iC(1:size(x_eig, 2));
        wh_to_delete = obj.xX.iC(size(x_eig, 2)+1:end);
        num_to_delete = length(wh_to_delete);
        
        obj.xX.X(:, wh_to_replace) = x_eig;
        obj.xX.X(:, wh_to_delete) = [];
        obj.xX.cond_assignments(wh_to_delete, :) = [];
        
        obj.xX.iB = obj.xX.iB - num_to_delete;
        obj.xX.iC = wh_to_replace;
        
end

end

