function subjects = canlab_list_subjects(basedir, varargin)
% subjects = canlab_list_subjects(basedir, wildcard1, wildcard2, etc.)
%
% Lists directories matching wildcard(s), returns cell array
% This cell array goes in EXPT.subjects in the CANlab standard
% variable/file structure.
%
% e.g.,
% subjects = canlab_list_subjects(pwd, 'nsf*');
%
% This lists directory names matching nsf* OR NSF*
% subjects = canlab_list_subjects(pwd, 'nsf*', 'NSF*');

d = [];
for i = 1:length(varargin)
    d = [d; dir(fullfile(basedir, varargin{i}))];
end

d(~cat(1, d.isdir)) = [];

d = char(d(:).name);
d = cellstr(d);
subjects = d';

end % function

