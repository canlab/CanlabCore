function data = montage_image_Worsley(image_name, varargin)
% ::
%
%    data = montage_image_Worsley(3D or 4D image name)
%
% Make a compact montage of some images
% Designed by Keith Worsley for pca_image.m
% Adapted by Tor Wager, Feb 2008
%
% Limits and Colormap are designed for component loadings between [-1 1]
%
% :Examples:
% ::
%
%    create_figure('Montage');
%    montage_image_Worsley('test_run1_pca.img');
%
%    data = montage_image_Worsley(imgname, 'pcacov')  % changes scaling and colormap
%    data = montage_image_Worsley(imgname, 'pcacov', [1 3 5])  % show only
%                                                           volumes 1, 3, 5 in image(s)

% ..
%    COPYRIGHT:   Copyright 2002 K.J. Worsley, 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that this copyright
%              notice appears in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
% ..
%
% ..
%    Tor Wager
% ..

scaletype = 'pcacor';  % Scale image and colors for PCA component viewing
whichvols = [];

for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'noscale', scaletype = 'none';
                    
                case 'pcacov', scaletype = 'pcacov';
                
                case 'whichvols', whichvols = varargin{i + 1};
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
end
    
% Get data
% -------------------------------------------------------------------------
image_name = expand_4d_filenames(image_name); % for SPM2 compatibility
V = spm_vol(image_name);

if isempty(whichvols)
    wh_vols = 1:length(V);
    
else
    % select only certain vols
    V = V(wh_vols);
end
    
p = length(V);
go_ok = 1;
if p > 10
    go_ok = input('More than 10 images! Are you sure you want to continue?'); 
end
if ~go_ok, return, end

data = spm_read_vols(V);


% Setup stuff
% -------------------------------------------------------------------------

numslices = V(1).dim(3);
% p = num comps
numys=V(1).dim(2);
numxs=V(1).dim(1);
   xr=1:numxs;
   yr=1:numys;
   %N=prod(V(1).dim(1:3));

nrow=round(sqrt(numslices/3.25/p));
nrow=max(nrow,1);

ncol=ceil(numslices/nrow);

numxr=length(xr);
numyr=length(yr);
bigmat=zeros(p*nrow*numyr,ncol*numxr);

% Create big mat
% -------------------------------------------------------------------------
r=0;
for k=1:p
   c=0; 
   for i=1:numslices
      bigmat((1:numyr)+r*numyr,(1:numxr)+c*numxr)=flipud(data(xr,yr,i,k)');
      c=c+1;
      if c>=ncol && i<numslices
         r=r+1;
         c=0;
      end
   end
   r=r+1;
end
xtick=((1:ncol*numxr)-1)/(ncol*numxr-1)*ncol-0.5;
ytick=((1:p*nrow*numyr)-1)/(p*nrow*numyr-1)*p+0.5;

switch scaletype
    case 'pcacor'
        zlim=[-1 1]; % [min(bigmat(:)) max(bigmat(:))];
    case 'pcacov'
        sdbm = nanstd(data(:));
        if(sdbm == 0)
            warning('Std dev of 0!');
            zlim = [-3 3];
        else
            zlim = [-3*sdbm 3*sdbm];
        end

    otherwise
        zlim = [min(bigmat(:)) max(bigmat(:))];
end

%
% Make image
% -------------------------------------------------------------------------
imagesc(xtick,ytick,bigmat,zlim);
xlabel('Slice (0 based)'); ylabel('Component');
switch scaletype
    case {'pcacor', 'pcacov'}
        title('Spatial components');
    otherwise
end
%colormap(spectral); colorbar;

set(gca,'YDir','Reverse')

% tor colormap stuff
% -------------------------------------------------------------------------
switch scaletype
    case {'pcacor'}
    cmap = colormap_tor([0 0 0], [1 1 0], [0 0 1], [0 1 1], [.5 .5 .5], [1 .4 .4], [1 0 0]);
    case {'pcacov'}
       cmap = colormap_tor([0 0 .5], [1 1 0], [0 0 1], [0 1 1], [0 0 0], [1 .4 .4], [1 0 0]);
     
    otherwise
     cmap = colormap_tor([0 0 0], [1 1 0], [0 0 1], [0 1 1], [.5 .5 .5], [1 .4 .4], [1 0 0]);
end       
colormap(cmap)
colorbar('vert')

axis tight

end
