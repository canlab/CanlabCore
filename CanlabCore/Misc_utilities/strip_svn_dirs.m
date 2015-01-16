% Removes the .svn dirs from the path
function strip_svn_dirs()
    strip_path_dirs('\.svn');
end