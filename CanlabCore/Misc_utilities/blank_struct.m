function C = blank_struct(A, blank_val)
% Utility function to create a blank structure, using an existing one
% as a template. 
%
% :Inputs:
%
%   **A:**
%        template structure
%
%   **blank_val:**
%        all fields will be set to this value, e.g. [] or ''
%
% :Output:
%
%   **C:**
%        structure containing fields of A , with values set to empty array
%
% Only works for single-element non-nested structures, for now
%
% ..
%    2/15/10 Joe Wielgosz
%    TODO: if B is array then loop..
% ..

fields = fieldnames(A);
C = struct;
for i = 1:length(fields)
    f = fields{i};
    C = setfield(C, f, blank_val);
    
end
