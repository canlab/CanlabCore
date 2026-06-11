function cl = region2struct(cl)
% region2struct Convert a region object to a simple struct array.
%
% Convert a region object to a simple MATLAB struct array, primarily for
% compatibility with other, older CANlab tools that expect a 'clusters'
% structure rather than a region-class object.
%
% :Usage:
% ::
%
%     cl = region2struct(cl)
%
% :Inputs:
%
%   **cl:**
%        A region-class object array.
%
% :Outputs:
%
%   **cl:**
%        A struct array containing the same fields as the input region
%        object, ready for use with legacy tools.
%
% :See also:
%   - cluster2region (the reverse transformation)
%   - region

warning off
for i = 1:length(cl)
    cl2(i) = struct(cl(i));
end
warning on

cl = cl2;

end

