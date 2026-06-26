function obj = removeblobs(obj)
% removeblobs Remove all rendered blobs from an fmridisplay object.
%
% Deletes blob graphics handles (and any associated legend handles) from
% every entry in obj.activation_maps, removes those entries from the
% object, and erases blobs from any registered surface objects via
% addbrain('eraseblobs', ...).
%
% :Usage:
% ::
%
%     obj = removeblobs(obj)
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object with one or more activation_maps and/or
%        surfaces attached.
%
% :Outputs:
%
%   **obj:**
%        The same fmridisplay object with deleted graphics handles and
%        with cleared entries removed from .activation_maps. Surfaces
%        in .surface are kept but have their blob colorings erased.
%
% :See also:
%   - fmridisplay
%   - addblobs
%   - removepoints
%   - addbrain

% Drop any views whose figures were closed, so we never try to erase blobs
% from deleted graphics handles.
obj = prune_dead_views(obj);

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

        myhan = myhan(ishandle(myhan));   % skip any handles whose figure was closed
        if isempty(myhan), continue, end

        myhan = addbrain('eraseblobs', myhan);

    end

end

% Keep an open controller panel in sync after layers are removed.
obj = update_controller(obj);

end
