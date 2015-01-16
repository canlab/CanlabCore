function history(dat)
% Display history for image_vector object

u = '_________________________________________________';
disp(u)
disp('History for image data object')
disp(u)

fprintf('Data description:\n%s\n', dat.dat_descrip);

disp(char(dat.history))

disp(u)


end
