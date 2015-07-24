function cl = atlas_spmanatomy_get_cluster(region_name)
% cl = atlas_spmanatomy_get_cluster(region_name)
%
% Returns a REGION object, cl, so that a named area from the SPM anatomy
% toolbox can be visualized.  
% 
% Requires SPM ANATOMY TOOLBOX V1.7 by Eickhoff (see code to update this.)
% 
% For a list of valid region names, do not enter region_name:
% atlas_spmanatomy_get_cluster
%
% Example:
% ba44_cl = atlas_spmanatomy_get_cluster('Area 44');
% 
% Visualize:
% cluster_orthviews(ba44_cl, {[0 0 1]}, 'add');
%
% Tor wager, Jan 2011

cl = region;

basename = 'AllAreas_v17_MPM.mat';
filename = which(basename);

if ~exist(filename, 'file')
    fprintf('Looking for atlas filename %s\nCannot find it.\n', basename)
    error('Exiting.')
end

load(filename)

imagename = which('AllAreas_v17.img');

if ~exist(imagename, 'file')
    fprintf('Looking for %s\nCannot find it.\n', imagename)
    error('Exiting.')
end

volInfo = iimg_read_img(imagename, 2);

names = char(MAP.name);


if nargin < 1 || isempty(region_name)
   disp('Run this function and enter a region name:') 
   disp(names);
   return
   
end

%% Get info about chosen region

wh = strmatch(region_name, names, 'exact');

if isempty(wh)
    disp('Warning! No match.')
    return
end

info = MAP(wh);


% Transfer info
N = fieldnames(cl);
M = fieldnames(info);
for i = 1:length(M)
    iname = M{i};
    if any(strmatch(iname, N))
        cl.(iname) = info.(iname);
    end
end

% Add other info

cl.M = volInfo.mat;
cl.dim = volInfo.dim;
cl.voxSize = diag(volInfo.mat(1:3, 1:3));
cl.center = mean(cl.XYZ', 1);
cl.mm_center = mean(cl.XYZmm', 1);
cl.source_images = imagename;
cl.numVox = size(cl.XYZ, 2);
cl.shorttitle = sprintf('%s', region_name);
cl.title = sprintf('%s from %s', region_name, filename);
cl.descrip1 = sprintf('%s, entry %3.0f in %s, mapped using %s', region_name, wh, filename, imagename);

end


