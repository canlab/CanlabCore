function m = celldat2matrix(dat)
% function m = celldat2matrix(dat)
%
% reshapes cell data in Chris format to a time x regions x subjects matrix
% concatenates across task states (conditions)
%
% cell data is {subjects x regions}, cells contain (time x condition)
% tor wager

for s = 1:size(dat,1)
    
    for r = 1:size(dat,2)
    
        x = dat{s,r}';
        x = x(:);
        
        m(:,r,s) = x;
        
    end
end

return


    
    
