function [poscl, negcl] = table(cl, varargin)
% Print a table of all regions in a region object (cl)
% 
% [poscl, negcl] = table(cl, [optional inputs])
% 
% Optional inputs:
% 'k'           : Print only regions with k or more contiguous voxels
% 'nosep'       : do not separate cl with pos and neg effects based on peak in .val
% 'names'       : name clusters before printing to table and output; saves in .shorttitle field
% 'forcenames'  : force naming of cl by removing existing names in .shorttitle field
%
% Outputs: 
% Returns region objects for cl with pos and neg effects, limited by size if entered
% and named if entered as optional input
%
% Copyright 2011, tor wager

k = 0;
dosep = 1;   
donames = 0; % name clusters before printing to table and output; saves in .shorttitle field
forcenames = 0; % force naming of cl by removing existing names in .shorttitle field

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case {'k', 'maxsize'}, k = varargin{i+1};
            case 'nosep', dosep = 0;
            case {'names', 'name', 'donames'}, donames = 1;
            case 'forcenames', forcenames = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


if k
    cl(cat(1, cl.numVox) < k) = [];
end

if donames
    
    if forcenames
        for i = 1:length(cl)
            cl(i).shorttitle = [];
        end
    end
    
    cl = cluster_names(cl);
end

if dosep
    % separate pos and neg
    [poscl, negcl] = posneg_separate(cl);
    fprintf('Positive Effects\n')
else
    % just return cl in poscl
    poscl = cl;
    negcl = [];
    fprintf('Table of all regions\n')
end

if ~isempty(poscl)
    cluster_table(poscl, 0, 0);
else
    disp('No regions to display');
end
fprintf('\nNegative Effects\n')

if ~isempty(negcl)
    cluster_table(negcl, 0, 0);
else
    disp('No regions to display');
end

end

