function Pso = getfullpath(Pspm)
% Searches for file in curr dir, then in specified dir
%
% :Usage:
% ::
%
%     Pso = getfullpath(Pspm)
%
% :Input:
%
%   **Pspm:**
%        is file name with no path or relative path
%
% :Output:
%
%   Returns absolute path name (full path)
%
% ..
%    tor wager
% ..

for i = 1:size(Pspm,1)
    
[d f e] = fileparts(deblank(Pspm(i,:)));
Ps = which(['./' f e]);

if isempty(Ps)
    
    if ~isempty(d), 
        cwd = pwd;
        try
            eval(['cd(''' d ''')']);,
        catch
            warning('Specified directory does not exist!')
            Ps = [];
            return
        end
        
        Ps = which(deblank(Pspm(i,:)));
        cd(cwd)
    end
end

if i > 1
    Pso = str2mat(Pso,Ps);
else
    Pso = Ps;
end

end

return
