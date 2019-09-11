function [str,sigmat] = correlation_to_text(a, cm, varargin)
% Print a matrix of correlation values, marked with * where significant
%
% :Usage:
% ::
%
%[str,sigmat] = correlation_to_text(matrix, value to mark with *,[cell array of names])
%
% input correlation matrix and critical value
% or raw data and empty cm value
% (determines threshold based on n)
%
% :Inputs:
%
%   **a:**
%        (1) A correlation matrix OR
%        (2) A matrix of data values to intercorrelate
%
%   **cm:**
%        If (1) above, a critical r value to mark with * in the output
%        If (2) above, empty (will calculate threshold)
%        OR
%        A matrix of significant values to mark with *
%
% see also print_correlation.m (deprecated), plot_correlation_matrix.m


if length(varargin) > 0
    
    nms = varargin{1}; 
    
else
    nms = {};
end

if ~iscell(nms)
    error('Names input must be a cell array.');
end

if ismatrix(cm)
    % We have entered a matrix to mark
    
    validateattributes(cm, {'numeric'}, {'square' '2d'});
    cm = logical(cm);
    
    str = sprintf('   \t');
    
elseif isempty(cm)

    % choose default value

    [rci,sig,z,p,cm] = r2z(.5,size(a,1), .05);

    a = corrcoef(a);

    str = sprintf('Crit=%3.2f\t', cm);

    cm = a > cm;  % logical matrix of which to mark
end

    

% names

for j=1:size(a,2)

    if length(nms) < j, nms{j} = ['V' num2str(j)]; end
    
    str=[str sprintf('%s\t',nms{j})];

end

str=[str sprintf('\n')];


% generate table

sigmat= false(size(a));

for i = 1:size(a,1)

    str=[str sprintf('%s\t',nms{i})];

    for j=1:size(a,2)

        if cm(i, j) && i ~= j, t='*';

            sigmat(i,j)=1;

        else,t='';
            
        end

        str=[str sprintf('%3.3f%s\t',a(i,j),t)];

    end

    str=[str sprintf('\n')];

end

disp(str)

end % function

