function [version_num, version_name, version_date] = check_matlab_version()
% [version_num, version_name, version_date] = check_matlab_version()
%
% version_date will be in datenum format
% Compare to a date like this: version_date > datenum('01-Jan-2016')
%
% [version_num, version_name, version_date] = check_matlab_version;
% if version_date > datenum('01-Jan-2016'), disp('You are running matlab 2016 or later!'), else, disp('old version!'), end

which_ver = ver('matlab');

version_num = str2double(which_ver.Version);

version_name = which_ver.Release;

version_date = datenum(which_ver.Date);


end

