%---------------------------------------------------------------------
% Utility function
% 2/15/10 Joe Wielgosz
%
% Create a blank structure, using an existing one as a template. 
%
% Outputs:
%   C: structure containing fields of A , with values set to empty array
%
% Inputs:
%   A: template structure
%   blank_val: all fields will be set to this value, e.g. [] or ''
%
% Only works for single-element non-nested structures, for now
function C = blank_struct(A, blank_val)
% TODO: if B is array then loop..
fields = fieldnames(A);
C = struct;
for i = 1:length(fields)
    f = fields{i};
    C = setfield(C, f, blank_val);
    
end