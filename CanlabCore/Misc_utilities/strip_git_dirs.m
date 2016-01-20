% Removes the .git dirs from the path
function strip_git_dirs()
    strip_path_dirs('\.git');
end
