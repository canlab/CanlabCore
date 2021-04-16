function obj = removeblobs(obj)

to_remove = [];

for i = 1:length(obj.activation_maps)
    
    if isfield(obj.activation_maps{i}, 'blobhandles')
        
        wh = ishandle(obj.activation_maps{i}.blobhandles);
        
        if any(wh)
            delete(obj.activation_maps{i}.blobhandles(wh))
            
            to_remove(end+1) = i;
        end
        
        if isfield(obj.activation_maps{i}, 'legendhandle')
            
            wh = ishandle(obj.activation_maps{i}.legendhandle);
            
            if any(wh)
                delete(obj.activation_maps{i}.legendhandle(wh))
            end
            
        end
    end
    
end

% Surfaces

obj.activation_maps(to_remove) = [];

if ~isempty(obj.surface)

    for i = 1:length(obj.surface)
        
        myhan = obj.surface{i}.object_handle;
        
        myhan = addbrain('eraseblobs', myhan);

    end
    
end


end
