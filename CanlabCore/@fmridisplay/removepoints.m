function obj = removepoints(obj)
% removepoints Delete plotted point handles from all montages in an fmridisplay.
%
% Iterates over each montage in obj.montage and deletes any valid
% graphics handles stored in the .plotted_point_handles field
% (created by addpoints).
%
% :Usage:
% ::
%
%     obj = removepoints(obj)
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object with points previously added via
%        addpoints.
%
% :Outputs:
%
%   **obj:**
%        The same fmridisplay object with .plotted_point_handles
%        deleted from each montage.
%
% :See also:
%   - fmridisplay
%   - addpoints
%   - removeblobs


for i = 1:length(obj.montage)
    
    if isfield(obj.montage{i}, 'plotted_point_handles')
        
        wh = ishandle(obj.montage{i}.plotted_point_handles);
        
        if any(wh)
            delete(obj.montage{i}.plotted_point_handles(wh))
        end
        
    end
    
end



end

