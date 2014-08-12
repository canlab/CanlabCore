function obj = removepoints(obj)


for i = 1:length(obj.montage)
    
    if isfield(obj.montage{i}, 'plotted_point_handles')
        
        wh = ishandle(obj.montage{i}.plotted_point_handles);
        
        if any(wh)
            delete(obj.montage{i}.plotted_point_handles(wh))
        end
        
    end
    
end



end

