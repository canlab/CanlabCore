function ok = checkMatlabVersion(datestr)
% :Usage:
% ::
%
%     ok = checkMatlabVersion(datestr)
% 
% check Matlab version against release date of version needed (or
% recommended)
%
% :Inputs:
%
%   **datestr:**
%        format is one of the forms taken by datenum.
%
% :Example: check for Matlab release date on or after Aug. 3, 2006
% ::
%
%    ok = checkMatlabVersion('8/3/2006') 
%
% ..
%    tor wager, 12/10/2006
% ..

    ok = 1;

% check Matlab version
        datethiscopy = datenum(version('-date'));
        dateneeded = datenum('8/3/2006'); % R2006a
        if datethiscopy < dateneeded
            disp('Warning: this version of Matlab is older than the recommended version.');
            disp('This function may not work correctly');
            ok = 0;
        end
        
end
