function [equality center] = doquality(Xcx,X)
% :Usage:
% ::
%
%     [equality] = doquality(Xcx,X)
%
% :Inputs:
%
%   **Xcx:**
%        binary indicator matrix of cluster assignments, 
%
%   **X:**
%        stimulus coordinates in group space
%
%        also: takes group spaces with zeros;


if size(Xcx>1);    
    for i = 1:size(Xcx,2)    % for each class
        tmp = mean(X(find(Xcx(:,i)),:),1);     % get center of this class
        center(i,:) = tmp;    
        % dist of all points in X from each center
        % number of cols of X is the number of dimensions; sums squared
        % vals across dims, takes sqrt to get Euclidean distance in N-d
        % space
        d(:,i) = sum((repmat(tmp,size(X,1),1) - X).^2,2).^0.5;
    end
    
    % rows of d are objects, columns are classes, values are dist from
    % class center
    for i=1:size(d,1)                   % i is the object (point)
        
        myclass = find(Xcx(i,:)==1);    % index of which class it is
        
        edist(i) = d(i,myclass);        %distance to center of own class
        
        otherclass = find(Xcx(i,:)==0);     % indices of columns for other classes
        
        otherdist(i) = min(d(i,otherclass));  %distance to center of cluster
        
    end

    % for each point, quality is distance to nearest neighbor - dist to own
    % cl / mak of those two
    equality = (otherdist - edist) ./ max([edist;otherdist]);
    
else
end
