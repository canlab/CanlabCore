function history(dat)
% history Display history for image_vector object.
%
% Prints the .dat_descrip and .history fields of an image_vector
% (or subclass) object to the command window, separated by underline
% rows.
%
% :Usage:
% ::
%
%     history(dat)
%
% :Inputs:
%
%   **dat:**
%        An image_vector (or subclass) object whose .history is a cell
%        or character array of provenance strings.
%
% :Outputs:
%
%   None. Output is printed to the command window.
%
% :Examples:
% ::
%
%     history(dat);
%
% :See also:
%   - descriptives
%   - image_vector

u = '_________________________________________________';
disp(u)
disp('History for image data object')
disp(u)

if ischar(dat.dat_descrip), fprintf('Data description:\n%s\n', dat.dat_descrip); end

disp(char(dat.history))

disp(u)


end
