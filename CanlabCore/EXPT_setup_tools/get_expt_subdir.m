function EXPT = get_expt_subdir(EXPT)
% EXPT = get_expt_subdir(EXPT)
%
% Gets names of individual subject directories and stores them in a field
% called EXPT.subjects
%
% example to create a new EXPT structure:
% EXPT = get_expt_subdir([]);


EXPT.subjects = [];
wc = input('Enter wildcard (e.g., 0*) :','s');
    
d = dir(wc);
    
for i = 1:length(d), EXPT.subjects{i} = d(i).name; end

return


    