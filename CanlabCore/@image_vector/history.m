function history(dat)
% Display history for image_vector object

u = '_________________________________________________';
disp(u)
disp('History for image data object')
disp(u)

if ischar(dat.dat_descrip), fprintf('Data description:\n%s\n', dat.dat_descrip); end

disp(char(dat.history))

disp(u)


end
