function clusters = sphere_roi_tool(varargin)
% :Usage:
% ::
%
%    clusters = sphere_roi_tool('mask','mask.img','bilat',1)
%    cl = sphere_roi_tool('mask',EXPT.mask,'bilat',0,'overlay',EXPT.overlay);
%
% A tool for building a number of spherical ROIs masked with some other
% mask (e.g., gray matter).  High-level end-user function.
%
% :Inputs:
%
%   arguments are in string - value pairs, in any order
%
%   **mask:**
%        mask image for spheres, default:
%        which('scalped_avg152T1_graymatter.img');
%
%   **bilat:**
%        make all rois bilateral, default is 1
%
%        radius of sphere is defined on a region-by-region basis
%
%   **write:**
%        followed by name of image to write out
%
% afterwards, try eliminating few-voxel regions:
% ::
%
%    cl(find(cat(1,cl.numVox)<10)) = [];
%
% and naming the clusters:
% ::
%
%    clusters = cluster_names(clusters);
%
% draw on existing regions:
% ::
%
%    cl = sphere_roi_tool('mask',EXPT.mask,'bilat',0,'overlay',EXPT.overlay,'add');
%
% use existing clusters
% ::
%
%    cl = sphere_roi_tool(cl,'mask',mask,'bilat',0,'overlay',ovl,'add');

% ..
%    defaults
% ..

bilat = 1;      % bilateral ROIs
rad = 8;        % radius in mm
donew = 1;      % make new orthviews
mask = which('scalped_avg152T1_graymatter.img');   % mask
overlay = which('scalped_avg152T1_graymatter.img');   % overlay
writename = 0;
includeall = 0;
includeexisting = 0;

% inputs
for i = 1:length(varargin)
    if strcmp(varargin{i},'mask'), mask = varargin{i+1};end
    if strcmp(varargin{i},'rad'), rad = varargin{i+1};end      % not used!
    if strcmp(varargin{i},'bilat'), bilat = varargin{i+1};end
    if strcmp(varargin{i},'overlay'), overlay = varargin{i+1};end
    if strcmp(varargin{i},'write'), writename = varargin{i+1};end
    if strcmp(varargin{i},'add'), donew = 0;end
    if isstruct(varargin{i}), existingcl = varargin{i}; includeexisting = 1;end
end

if isempty(overlay), warning('cannot find overlay.'); overlay = which('scalped_avg152T1_graymatter.img'); end
if isempty(overlay),error('Cannot find either EXPT.overlay or default overlay.'); end

% initial views
try
    %cl(1).mm_center = [0 0 0]; cl(1).XYZmm = [0 0 0]'; cl(1).XYZ = cl(1).XYZmm;cl(1).Z = 1;cl(1).title = '';
    %icbm_localize(cl)

    maskcl = mask2clusters(mask);
    % try to make clusters faint
    for i = 1:length(maskcl)
        maskcl(i).Z = .1 * ones(1,size(maskcl(i).XYZ,2));
        maskcl(i).Z(1) = 1;
    end
    if donew, cluster_orthviews(maskcl,{[1 0 0]},'overlay',overlay); end
    %spm_image('init',overlay);
    spm_orthviews('Reposition',[0 0 0])

catch
    % should give an error, as we haven't fully defined clusters
end

if includeexisting
   % inmask = clusters2mask(existingcl,dims);  %spm_read_vols(spm_vol(mask));
   disp('Found existing clusters.')
    clusters = existingcl;
else
    clusters = [];
end
    
go = 1;
r = 10;

rem_clusters = [];
[volInfo,maskindx] = iimg_read_img(mask,1);

% set up callback that gets cluster when you click
data.volInfo = volInfo;
data.r = r;
han = findobj('Tag','Graphics');
if isempty(han), error('SPM orthviews did not open for some reason?'); end
guidata(han,data)
%fhan = @get_cl;
set(han,'WindowButtonUpFcn','try, cl = get_cl; catch, disp(''Press and hold a sec to get cluster.''), end');
rehash

% instructions
fprintf(1,'Enter radius in mm to set a new radius, 1 to save the current cluster, or 0 to exit.\n');

while go

    val = input('0 to exit / 1 to add / 999 to remove / or enter new radius? ');

    if isempty(val)
        % do nothing
    else
        switch val
            case 0, go = 0;

            case 1
                fprintf(1,'\nSaving the cluster listed above.\n');
                data = guidata(han);
                cluster_orthviews(data.cl,{[1 1 0]},'add');
                try, clusters = merge_clusters(clusters,data.cl); catch, disp('Error adding cluster. Empty?'), end

                
            case 999
                fprintf(1,'\nRemoving the cluster listed above.\n');
                data = guidata(han);
                cluster_orthviews(data.cl,{[0 0 .2]},'add');
                try, rem_clusters = merge_clusters(rem_clusters,data.cl); catch, disp('Error adding cluster. Empty?'), end
                
            otherwise
                data = guidata(han);
                data.r = val;
                guidata(han,data)
                fprintf(1,'\nRadius is now %3.0f mm\n',data.r);
                

        end
    end
end

if bilat
    %pos2 = pos; pos2(1) = -pos2(1); pos = [pos; pos2];    % reverse x coord to make bilateral (bilat mode)
end

fprintf(1,'Finishing up: combining included and excluded.\n');
dims = volInfo.dim(1:3);
exmask = zeros(dims);

% concatenate clusters into contig regions
%[tmp,tmp,clusters] = clusters2mask(clusters,clusters(1),0,mask);
if ~isempty(clusters)
    inmask = clusters2mask(clusters,dims);
    
    if writename
        V = volInfo; V.fname = 'include.img';
        spm_write_vol(V,inmask);
        disp(['Written: ' V.fname])
    end
end

if ~isempty(rem_clusters)
    exmask = clusters2mask(rem_clusters,dims);
    
    if writename
        V = volInfo; V.fname = 'exclude.img';
        spm_write_vol(V,exmask);
        disp(['Written: ' V.fname])
    end
end

if exist('inmask', 'var') && exist('exmask', 'var')
    finalimg = inmask & ~exmask;
elseif exist('exmask', 'var')
    finalimg = ~exmask;
elseif exist('inmask', 'var')
    finalimg = inmask;
else
    error('???')
end

clusters = mask2clusters(finalimg,volInfo.mat);

if writename
    V = volInfo; V.fname = dowrite;
    spm_write_vol(V,finalimg);
    disp(['Written: ' V.fname])
end
    
set(gcf,'WindowButtonUpFcn','');


return

% subfunction: Window button up

function cl = get_cl

% get data
han = findobj('Tag','Graphics');
if isempty(han), error('SPM orthviews did not open for some reason?'); end
data = guidata(han);
volInfo = data.volInfo;
r = data.r;

% get coordinates (mm and voxel)
pos = round(spm_orthviews('Pos')');

xyz = mm2voxel(pos,volInfo);

% voxel radius
rv = round(mean(abs(r ./ diag(volInfo.mat(1:3,1:3))')));   % voxel coords

indx = iimg_xyz2spheres(xyz,volInfo.xyzlist,rv);

fprintf(1,'\nCenter: %3.0f,%3.0f,%3.0f.  Radius: %3.0f mm.  Voxels: %3.0f',pos(1),pos(2),pos(3),r,sum(indx));

cl = iimg_indx2clusters(indx,volInfo);

% add to gui data because i can't seem to get it to return
data.cl = cl;
guidata(han,data);

return
