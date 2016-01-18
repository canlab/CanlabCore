function strip_git_dirs()
% Removes the .svn dirs from the path
    strip_path_dirs('\.git');
end
