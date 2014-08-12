function rename_lowercase
% rename_lowercase
%
% tor wager, renames all files in dir with lower case versions
% if all characters are capitalized.
%
% batch example:
%d = dir('0*'); for i=1:length(d),cd(d(i).name),rename_lowercase, cd .., end

d = dir; for i = 3:length(d), 
    nm1 =(d(i).name);nm2=upper(nm1);
    str=['!mv ' nm1 ' ' nm2]; 
    
    if strcmp(nm1,'SPM.MAT'),
        str = ['!mv ' nm1 ' SPM.mat'];
        disp(str)
        eval(str)
    elseif strcmp(lower(nm1),nm1)
        disp(str)
        eval(str)
    end
end

return