function x_cell = canlab_mat2cell(x, s)
% Simpler replacement for mat2cell for common use cases in our tools.
%
% :Usage:
% ::
%
% x_cell = canlab_mat2cell(x, s)
%
% :Inputs:
%
%   **x:**
%        m x k matrix to be divided into cells
%
%   **s:**
%        vector of integers 1...n describing which rows of x belong to each cell
%        rows with 1 go in cell 1, 2 in cell 2, etc. 
%        missing integers will result in empty cells
%        rows with 0 will not be placed into cells
%        
% :Outputs:
%
%   **x_cell:**
%        cell array, {1 x n}, with x submatrices
%
% :Examples:
% ::
%
% Define integer list of subjects s and obs x variables matrix s
% s = ST_cleaned.subj_idx;
% x = ST_cleaned.fine_regions;
% x_cell = canlab_mat2cell(x, s);
% 
% Now calculate correlation matrices for each subject
% out = cellfun(@corr, x_cell, 'UniformOutput', false);
%
% Spearman's correlation instead:
% out = cellfun(@(x) corr(x, 'Type', 'Spearman'), x_cell, 'UniformOutput', false);
%
% Concatenate for further testing:
% out = cat(3, out{:});

%u = unique(s);     % this would prevent missing integers from resulting in empty cells

% Validate input attributes
validateattributes(s,{'numeric'},{'2d' 'integer' 'column'});
validateattributes(x,{'numeric'},{'2d' 'finite'});

u = 1:max(s);

n = length(u);

x_cell = cell(1, n);

for i = 1:n
    
    wh = s == u(i);
    
    x_cell{i} = x(wh, :);
    
end

end % function
