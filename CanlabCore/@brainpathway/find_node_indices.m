
function to_extract = find_node_indices(obj, varargin)
% Identify nodes to extract from brainpathway object given an integer
% vector of which nodes to select, cell array of strings for node names,
% or combination of the two. Returns a logical vector of nodes in set
% 
% obj : a brainpathway object
% 
% Optional inputs:
% a cell array of strings containing node names (matched to obj.node_labels)
% an integer vector of which nodes
%
% 'flatten' (not used here - would collapse atlas regions if used in select_atlas_regions)
%
% Outputs:
% to_extract: a logical

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

strings_to_find = [];
integers_to_find = [];
doflatten = false;

% optional inputs with default values
for i = 1:length(varargin)
    
    if iscell(varargin{i})
        strings_to_find = varargin{i};
        
    elseif isnumeric(varargin{i})
        integers_to_find = varargin{i};
        
    elseif ischar(varargin{i})
        switch varargin{i}
            
            case 'flatten', doflatten = true;
                
                %             case 'xxx', xxx = varargin{i+1}; varargin{i+1} = [];
                
            otherwise
                warning(['Unknown input string option:' varargin{i} '. Assuming it might be an atlas label. Place atlas labels in a cell array']);
                strings_to_find{1} = varargin{i};
        end
    end
end


% -------------------------------------------------------------------------
% INIT
% -------------------------------------------------------------------------

k = length(obj.region_indx_for_nodes); % size(obj.node_dat, 2);
to_extract = false(1, k);

% -------------------------------------------------------------------------
% FIND BY STRING
% -------------------------------------------------------------------------

for i = 1:length(strings_to_find)
    
    % Find which names match
    wh = ~cellfun(@isempty, strfind(obj.node_labels, strings_to_find{i}));
    
    to_extract = to_extract | wh;
    
end

% -------------------------------------------------------------------------
% FIND BY NUMBERS
% -------------------------------------------------------------------------

to_extract(integers_to_find) = true;

if ~any(to_extract)
    error('No nodes identified to extract.');
end

end % subfunction
