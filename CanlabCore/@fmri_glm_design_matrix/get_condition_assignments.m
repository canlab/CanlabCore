function cmatrix = get_condition_assignments(obj)
% Condition assignments
%   - Indicator matrix coding for which columns in X belong to the same
%     modeled condition, and are part of the same HRF fit
%   - There is one set of columns for each condition modeled, and one set of
%     columns for each parametric modulator of each condition
%   - Because parametric modulators may not exist for all conditions, we need
%     to build this dynamically for modulators.
%
% Design matrix build (which calls method get_session_X) builds columns in
% this order:
%
% All within Session:
% Regressors of interest, basis functions within conditions
% Parametric modulators, basis functions within conditions
% Covariates of no interest
%
% Then:
% Baselines (session/run intercepts)
%
% This method is called automatically in the build method.

switch obj.build_method
    case 'Separate sessions'
        nsess = length(obj.Sess);


    case 'Concatenated sessions'
        nsess = 1;

        
    otherwise
        error('Unknown build method for fmri_model object. See help for valid strings.');
        
end


for s = 1:nsess
    nconds(s) = length(obj.Sess(s).U);
end
if any(diff(nconds)), error('Variable numbers of conditions; this should have been caught in build.'); end


% cell of indicator matrices for each session separately
all_indic = cell(1, nsess);  

for s = 1:nsess
    
    % regs of interest
    % -----------------------------------------------------------
    sessindic = {};
    
    for i = 1:nconds(s)
        
        
        % Allow for variable number of bfs.
        nbf = size(obj.xBF(i).bf, 2);

        % indicate a set of regressors based on basis set
        sessindic{i} = true(nbf, 1);
        
    end
    
    % concatenate into block diagonal
    sessindic = blkdiag(sessindic{:});
    
    % Modulators
    % -----------------------------------------------------------
    pmindic = {};
    
    for i = 1:nconds(s)
        is_pm = isfield(obj.Sess(s).U(i), 'P') && ~isempty(obj.Sess(s).U(i).P) && isfield(obj.Sess(s).U(i).P, 'P') && ~isempty(obj.Sess(s).U(i).P.P);
        
        if is_pm
            % Allow for variable number of bfs.
            nbf = size(obj.xBF(i).bf, 2);
        
            pmindic{end + 1} = true(nbf, 1);
        end
    end
    if ~isempty(pmindic)
        pmindic = blkdiag(pmindic{:});
    else
        pmindic = [];
    end
    
    % add them together in block diag form
    all_indic{s} = blkdiag(sessindic, pmindic);
    
end % session

cmatrix = logical(blkdiag(all_indic{:})); 

ncols = size(obj.xX.X, 2);

if size(cmatrix, 1) > ncols
    error('More conditions of interest than rows in obj.xX.X!! Must build xX before running this.');
end

% add rows of zeros to match total number of columns
z = false(ncols - size(cmatrix, 1), size(cmatrix, 2));

cmatrix = [cmatrix; z];


end % function
