function myLight = standardMRIlighting(option,handles)
% :Usage:
% ::
%
%    myLight = standardMRIlighting(option,handles)
%
% :Inputs:
%
%   **option:**
%        'full' - all lighting adjustments
%
%        'reflectance' - ambient strength and reflectance only
%
%   **handles:**
%        [isosurfaceHandle isocapsHandle]            

if strcmp(option,'full')
    
    view(135,30) 

    myLight = lightangle(45,30); 
    %set(gcf,'Renderer','zbuffer'); 
    lighting phong
    set(myLight,'Tag','myLight')
end

if ~isempty(handles)
    if length(handles) > 1
        try
            
        set(handles(2),'AmbientStrength',.6)
        catch
        end
    end
    
    try
        set(handles(1),'SpecularColorReflectance',0,'SpecularExponent',50)
    catch
    end
end

if strcmp(option,'full')
    % this is the key to avoiding dark surface head.
    [az,el] = view;
    myLight = lightangle(az,el);
else myLight = [];
end

axis image

drawnow

return
