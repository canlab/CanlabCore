%------------------------------------------------------------------------------
% Sub-function for imaging patch
%------------------------------------------------------------------------------
function [p1,p2] = imagePatch(X,Y,Z,D,Ds,inValues)

    if isempty(inValues)
        inValues = [nan nan nan nan nan nan];
    end
    
    if isempty(X) | isempty(Y) | isempty(Z)
        % no xyz coordinates
        FV2 = isosurface(Ds,50);
        IS2 = isocaps(D,50);

        % Draw figure patches
        p1 = patch(FV2,'FaceColor',[1,.75,.65],'EdgeColor','none');
        p2 = patch(IS2,'FaceColor','interp','EdgeColor','none');
    else
        
        % Define subvolume for cutaway view
        [x y z A] = subvolume(X,Y,Z,D,inValues);
        [x y z As] = subvolume(X,Y,Z,Ds,inValues);
        FV2 = isosurface(x,y,z,As,50);
        IS2 = isocaps(x,y,z,A,50);

        % Draw figure patches
        p1 = patch(FV2,'FaceColor',[1,.75,.65],'EdgeColor','none');
        p2 = patch(IS2,'FaceColor','interp','EdgeColor','none');
        % isonormals(A,p1)
    end
    
    drawnow
return