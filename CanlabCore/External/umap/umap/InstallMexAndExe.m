%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function InstallMexAndExe
[mexFileName, mexFile, umapFolder]=UmapUtil.LocateMex('sgd');
if exist(mexFile, 'file')
    warning(['You already have the mex file: ' mexFile]);
    fprintf('If you wish to rebuild then first remove this!\n\n');
    return;
end
disp('This allows you to invoke run_umap with ''method''=''MEX'' which');
disp('provides our fastest version of stochastic gradient descent');
fprintf('... otherwise it is done slower with Java\n\n');
cppFile=fullfile(umapFolder, 'sgdCpp_files', 'mexStochasticGradientDescent.cpp');
if ~exist(cppFile, 'file')
    error('Can not find the C++ file ... sigh');
end
try
    mex(cppFile)
catch ex
    ex.getReport
    msg('Setup C++ compiler correctly with "mex -setup"!');
end

if ~exist(mexFile, 'file')
    m=fullfile(pwd, mexFileName);
    if exist(m, 'file')
        movefile(m, umapFolder);
    end
end
if exist(mexFile, 'file')
    fprintf(['\nGOOD...The build created "' mexFileName ...
        '"\n\t in folder "' umapFolder '"\n' ]);
    fprintf('\nIf you distribute this to other computers, remember to sign it!!\n');
    fprintf('For example,at the Herzenberg lab on MAC computers we type:\n');
    fprintf('codesign -s "Stanford University" mexStochasticGradientDescent.mexmaci64\n\n');
else
    warning(['Could not build ' mexFile '???']);
end
disp('The exe is the same C++ for doing stochastic gradient descent. But ');
disp('it runs slower. We keep it available for MatLab accleration education');
disp('so that you can build it with CLang++ and the build script in the ');
fprintf('sgdCpp_files sub folder of umap... then invoke run_umap ''method''=''C++'' !\n\n');
fprintf('Have fun reducing with UMAP!!\n\n');
end