%---------------------------------------------------------------------
% Utility function
% 2/15/10 Joe Wielgosz
%
% Combine the fields of two structures. 
%
% Outputs:
%   C: structure containing fields of both A and B
%
% Inputs:
%   A, B: structs to be combined
%   prefix: string applied to fields of structure B
%
% Only works for single-element non-nested structures, for now
function C = combine_structs(A, B, prefix)
C = A;

% TODO: if B is array then loop..
fields = fieldnames(B);
for i = 1:length(fields)
    if exist('prefix', 'var')
        f = [prefix fields{i}];
    else
        f = fields{i};
    end
    val = getfield(B, fields{i});
    % TODO: check for conflicts...
    % TODO: if substruct then recursively combine subfields..
    C = setfield(C, f, val);
    
end