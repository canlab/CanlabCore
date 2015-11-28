function cl = cluster2subclusters(cl_in,class)
% Take a single cluster cl_in and separate into subclusters based on
% vector of integers class
%
% Class must code unique subclusters subcluster order is only preserved
% if class contains all integers from 1 to nclasses:
%
% i.e., class 3 will only be in subcluster 3 if there are no missing class
% numbers in class
%
% ..
%    tor wager, july 06
% ..

classes = unique(class);    % vector of class numbers

nclust = length(classes);

% Fill in fields
for i = 1:nclust

    % voxels in this cluster
    wh_incluster = find(class == classes(i));
    nvox = length(wh_incluster);

    cl(i).title = sprintf('Subcluster of %3.0f vox from %s',nvox, cl_in.title);

    % copy these fields
    for f = {'threshold' 'M' 'dim' 'voxSize'}
        if isfield(cl_in,f{1}), cl(i).(f{1}) = cl_in.(f{1}); end
    end
    
    % add new values for these fields
    cl(i).name = '';
    cl(i).numVox = nvox;
    
    % select in_subcl voxels for these fields
    for f = {'XYZ' 'XYZmm' 'Z'}
        cl(i).(f{1}) = cl_in.(f{1})(:,wh_incluster);
    end
    
    % center of mass
    cl(i).mm_center = center_of_mass(cl(i).XYZmm,cl(i).Z);
    
end
