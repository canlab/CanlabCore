% just run it.
%
% Prompts you for all necessary info, including image to draw on
% and output image file name.
%
% Writes a mask *img file in analyze format that you can view in spm
%
% Functions called:
% c:\tor_scripts\voistatutility\readim2.m
%   read_hdr.m
%
% c:\tor_scripts\roiutility\imgmanip\getvoxels3.m
%   C:\matlabR12\toolbox\matlab\elmat\ind2sub.m
%
% c:\tor_scripts\roiutility\cluster\tor_extract_rois.m
%
% C:\matlabR12\toolbox\matlab\datatypes\squeeze.m

% c:\tor_scripts\voistatutility\nanmean.m
% center_of_mass.m

% C:\matlabR12\toolbox\matlab\elmat\repmat.m
% (calls other spm functions)
%
% c:\tor_scripts\fixedfxutility\voxel2mm.m
%
% spm_vol.m
% spm_read_vols.m
%
% ..
%    by Tor Wager  last modified 9/22/02
% ..
%
% :Note: readim2 works on little endian systems (Linux).
% If you're using Unix, you should modify this script
% to call the function readim2_b.m instead.
%

P = spm_select(1,'image','Choose whole-brain anatomical mask');
%P = spm_get(1,'*.img','Choose whole-brain anatomical mask');


V = spm_vol(P);
vol = spm_read_vols(V);
vol = double(vol);


dosagg = input('Press 1 to choose saggital slices, or 0 for axial.');
% if saggital, use vol, not rvol
% y is not reversed
% ginput x = brain y, ginput y = brainz, slice = x

if dosagg
    [rvol2] = readim2(vol,'p','sagg');
else
    for i = 1:size(vol,3)
        rvol(:,:,i) = rot90(vol(:,:,i));
    end

    [rvol2] = readim2(rvol,'p');
end

colormap gray
wslices = input('Enter range of slices (e.g, 20:30) ');
close

if dosagg
    [rvol2,hdr,h] = readim2(vol,'p','sagg','noflipy',wslices);
    colormap gray
    [voxels,mask] = getvoxels3(h,vol,wslices,'sagg');
else
    [rvol2,hdr,h] = readim2(rvol,'p','ax','noflipy',wslices);
    colormap gray
    [voxels,mask] = getvoxels3(h,vol,wslices);
end


%for i = 1:size(mask,3)
%    rmask(:,:,i) = rot90(mask(:,:,i),3);
%end

V.fname = input('Enter output filename, no quotes: ','s');
spm_write_vol(V,mask);
[d,fname,e] = fileparts(V.fname);
saveas(gcf,fname,'fig')

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% cluster structure stuff for ROIs
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CLU.M = V.mat;
CLU.voxSize = diag(V.mat)'; 
CLU.voxSize = CLU.voxSize(1:3);
CLU.VOX = CLU.voxSize;
CLU.XYZ = voxels';
CLU.XYZmm = voxel2mm(voxels',V.mat);
CLU.Z = ones(1,size(CLU.XYZ,2));
CLU.crit_t = 1;
CLU.cl_size = 0;
CLU.u = CLU.crit_t;
CLU.k = CLU.cl_size;
CLU.title = fname;

clusters = tor_extract_rois([],CLU,CLU);
eval(['save ' fname '_clusters CLU clusters'])

disp(' ')
newim = input('Press 1 to choose a display image, or anything else to use your prior mask.');
if newim, P = spm_get(1,'*.img','Choose overlay image');,end

spm_image('init',P)
for i = 1:length(clusters)
    spm_orthviews('AddColouredBlobs',1,clusters(i).XYZ,clusters(i).Z,CLU.M,rand(1,3))
end
