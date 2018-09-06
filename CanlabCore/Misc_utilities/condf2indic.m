function [indic, xlevels] = condf2indic(X, varargin)
% Create N x k indicator matrix of ones/zeros from N x 1 list of numeric
% values (X)
%
% :Usage:
% ::
%
%     [indic, xlevels] = condf2indic(X, ['integers'])
%
% :Inputs:
%
%   **'integers':**
%       Treat X as integer vector. Check that all are integers, and return
%       one column per integer in order, excluding zeros. Insert NaNs in
%       columns whose corresponding integers are missing.  Returns an n x k
%       matrix indic, where n is size(X, 1) and k is max(X).
%
%       'integers' may be followed by another variable specifying the
%       number of columns to use. this is useful because integers at the
%       end of a list in X may be missing; this will append columns for
%       those.
%
% :Outputs:
%
%   **indic:**
%        is returned as single precision type, because it can then be
%        used as a design matrix in a GLM
%
%   **xlevels:**
%        are the values of X corresponding to columns of indic
%
% :Examples:
% condf = [1 1 1 2 2 2 4 4 4 6 6 1 1 1 8 8]';
% [indic, xlevels] = condf2indic(condf); indic, xlevels
% [indic, xlevels] = condf2indic(condf, 'integers'); indic, xlevels
% [indic, xlevels] = condf2indic(condf, 'integers', 10); indic, xlevels

% ..
%    tor wager, nov 2007. Edited Feb 2018 to expand functionality to
%    'integers' case, for use with atlas object.
%    Aug 2018: fixed some inconsistencies in 'integers' option.
% ..

% PRELIM CALCULATIONS

X = double(X); % for objects/special types

[xlevels, ~, indx] = unique(X);

max_integers = length(xlevels);


% INPUTS

dointegers = false;

if any(strcmp(varargin, 'integers'))
    
    dointegers = true;
    
    wh_input = strcmp(varargin, 'integers');
    
    if length(varargin) > find(wh_input)      % we have entered max integers
        
        max_integers = varargin{wh_input + 1};
        
    else
        
        max_integers = max(xlevels);
        
    end
    
end

% INITIALIZE - cols equal to num unique entries, or all integers

indic = single(false(length(indx), max_integers));

if dointegers
    
    % Check integers
    u = unique(X);
    if ~all(u == round(u)) || any(u < 0)
        warning('Some condition function values are not integers.');
    end
    
    % Remove zero-valued integer
    wh = find(xlevels == 0);
    xlevels(wh) = [];
    indic(:, wh) = [];
    if ~isempty(wh)
        indx(wh) = 0; 
        indx(indx >= wh) = indx(indx >= wh) - 1; % if there was a zero, then subtract 1 from each indx value for conditions greater than the one dropped
    end
    %if ~isempty(wh), indx = indx - 1; end  % if there was a zero, then assume 0 is first indx value 
    
    for i = 1:length(xlevels)
        indic(indx == i, xlevels(i)) = 1;
    end
    
    xlevels = (1:max_integers)';

else
    % no integers, do not preserve empty columns
    
    for i = 1:length(xlevels)
        indic(indx == i, i) = 1;
    end

end % integers

end % main function
