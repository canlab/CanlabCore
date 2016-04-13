function [array,hdr,h,whichslices,rows,cols,figh] = readim2(varargin)
% :Usage:
% ::
%
%    [array,hdr,h,whichslices,rows,cols,figh] = readim2(basename or array [opt],'p' [opt], 'sagg' or 'cor' [opt],flipy[opy],range [opt])
%
% :Inputs:
%
%   - basename of file, without image extension ,OR
%   - 3-D array in the workspace to plot, OR
%   - nothing, to browse for file
%
%   **p:**
%        to plot montage of slices to the screen
%
%   **sagg:**
%        to rotate to saggital view
%
%   **cor:**
%        to rotate to coronal view
%
%   **t:**
%        to save array as double instead of int16 - to save negative t values.
%
% :Outputs:
%
%   - 3d array of image
%   - hdr of image
%   - handles for axes of montage if plotting
%
% :Special Features:
%
%   **range:**
%        specified in mm, must be LAST and FIFTH input argument. 
%          - OR range can specify slices, e.g. 1:4:28
%
%   **clim:**
%        color limits for axis plot, must be 1st or 2nd argument.  form: [-1 1]
%
% ..
%   Adapted 12/06/00 by Tor Wager from Luis Hernandez' original script
%   Last modified 10/19/01 by Tor
% ..

warning off 

if nargin < 1
	arg1 = [];arg2 = [];arg3 = [];arg4 = [];
else [array,basename,orient,print,tmap,flipy,usespm] = setup(nargin,varargin{1:end});
	justbase = basename;
end

if size(size(array),2) == 3
	basename = 'wkspace variable';
elseif nargin < 1 | ~exist('basename')
    	[filename,pathname] = uigetfile('*.img','pick an img file.');
    	basename = [pathname filename(1:end-4)];
	justbase = filename(1:end-4);
elseif isempty(basename) | strcmp(basename,'null') | ~ischar(basename)
    	[filename,pathname] = uigetfile('*.img','pick an img file.');
    	basename = [pathname filename(1:end-4)];
	justbase = filename(1:end-4);
end

% ----------------------------------------------------------------------------------
% * check to see if it exists, if not look somewhere else. this is totally unnecessary
% ----------------------------------------------------------------------------------
% first try taking off the .img extension.
if ~(exist([basename '.img']) == 2)
	basename = basename(1:end-4);
end
	
if ~(exist([basename '.img']) == 2)
	warning('File cannot be found in specified directory.')
	basename = which([justbase '.img']);
	basename = basename(1:end-4);
	
	disp(['Using: ' basename])
end


% LOADING THE ARRAY
if isempty(array)
	% Luis hernandez
	% last edit 6-29-2000
	% get the header for the basename
	eval(['hdr = read_hdr(''' basename '.hdr'');'])
    eval(['[pFile,messg] = fopen(''' basename '.img'',''r'');']) 
 	if pFile == -1
      		disp(messg);   
    		return;
   	end

	switch hdr.datatype     
   	case 0
      		fmt = 'int8';
   	case 2
      		fmt = 'uint8';
   	case 4
      		fmt = 'int16'; %'short';
   	case 8
      		fmt = 'int';
   	case 16
      		fmt = 'float';
   	case 32
      		fmt = 'float';
      		xdim = hdr.xdim * 2;
      		ydim = hdr.ydim * 2;
   	case 64
      		fmt = 'uint64';
   	otherwise
         	warning(['Data Type ' num2str(hdr.datatype) 'Unsupported.  Switching to big-endian...']);
	 	eval(['hdr = read_hdr_b(''' basename '.hdr'');'])
		status = fclose(pFile);
    		eval(['[pFile,messg] = fopen(''' basename '.img'',''r'',''b'');']) 
		if pFile == -1
      			disp(messg);   
    			return;
   		end

		switch hdr.datatype     
   		case 0
      			fmt = 'int8';
   		case 2
      			fmt = 'uint8';
   		case 4
      			fmt = 'int16'; %'short';
   		case 8
      			fmt = 'int';
   		case 16
      			fmt = 'float';
   		case 32
      			fmt = 'float';
      			xdim = hdr.xdim * 2;
      			ydim = hdr.ydim * 2;
   		case 64
      			fmt = 'uint64';
   		otherwise
			error('Cannot read native or big-endian format.')
		end
	end

    % ===== make space for double if tmap, int16 otherwise =====
    if tmap, array = zeros(hdr.xdim,hdr.ydim,1);
    else array = int16(zeros(hdr.xdim,hdr.ydim,1));
    end
    warning off
	for i=1:hdr.zdim    
      a  = (fread(pFile,[hdr.xdim, hdr.ydim], fmt)) ;
        a(isnan(a)) = 0;
		array(:,:,i) = a; 
        if mod(i,50) == 0,disp(['    loaded ' num2str(i) ' slices']),end
    end
    warning on
	fclose(pFile);
	
else		% if it's not a string
   hdr.zdim = size(array,3);
   justbase = 'workspace variable';
end

%ROTATIONS
if flipy
	disp('flipping y dimension of img (columns of array).')
	array = flipdim(array,2);
end

if strcmp(orient,'sagg')
	array = permute(array,[3 2 1]);
	array = flipdim(array,1);
elseif strcmp(orient,'cor')
	array = permute(array,[3 1 2]);
	array = flipdim(array,1);
end

% restrict range, if necessary
if nargin > 4 & ~ischar(varargin{5}) %  & sum(size(varargin{1})==[1 2]) == 2
   range = varargin{5};
   if size(range,2) == 2 % range in mm, convert to slices
   		hdr.effzdim = size(array,3);
      slicerange = round(range(1)/hdr.effzdim):round(range(2)/hdr.effzdim); 
   else	% range should be row vector of slices
      slicerange = range;
   end
      array = array(:,:,slicerange);
end

[rows,cols,skip,whichslices] = getplotdims(array,print);
	
h = [];
%DISPLAYING IMAGE
if print
   whos array 
   array2 = double(array);
   
   if usespm
      if(findobj('Tag','Graphics')),figure(findobj('Tag','Graphics'))
      else spm_figure('Create','Graphics','Graphics')
      end
      figh = gcf; clf;
   else
   		figh = figure;
   end
      
    zoom = 6/max(rows,cols);
    set(gcf,'Position',[46 16 1106 916]);
    if size(varargin{1},2) == 2 %sum(size(varargin{1})==[1 2]) == 2 
       clim = varargin{1};
    elseif size(varargin{2},2) == 2 %sum(size(varargin{2})==[1 2]) == 2 
       clim = varargin{2};
    else clim = [min(min(min(array2))) max(max(max(array2)))];
    end
	if clim(1) == clim(2), clim(2) = clim(2) + 1;,warning('All elements have same value!'),end
	index = 1;
	for i = whichslices	
      h(index) = subplot(rows,cols,index);
  		imagesc(array(:,:,i),[clim]);
  		axis off
  		axis square
  		if rows < 6 & cols < 6	 
		 	j = get(gca,'Position'); j = [j(1:2) .12 .12];
  		 	set(gca,'Position', j);
		end
		index = index + 1;
        camzoom(zoom)
        dasp = get(gca,'DataAspectRatio');
        switch orient
        case 'sagg'
            daspect([dasp(1)*.78 dasp(2) 1])
        end
	end
	h(index) = subplot(rows,cols,index);
   imagesc(clim); axis off
   colorbar('horiz') % set(H,'XTick',[clim(1) clim(2)])
	cla,
   text(0,1,justbase);
   text(0,1.3,['Range = ' num2str(min(min(min(array2)))) ' to ' num2str(max(max(max(array2))))])
end

warning on

return

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %

function [array,basename,orient,print,tmap,flipy,usespm] = setup(nargin,varargin)
   
print = 0;tmap=0;orient='ax';flipy=0;usespm=0;array=[];basename='null';

for i = 1:nargin
   if strcmp(varargin{i},'p'),print = 1;,end
   if strcmp(varargin{i},'t'),tmap = 1;,end
	if strcmp(varargin{i},'flipy'),flipy = 1;,end
   if strcmp(varargin{i},'spm'),usespm = 1;,end
   if strcmp(varargin{i},'cor'),orient = 'cor';,end
	if strcmp(varargin{i},'sagg'),orient = 'sagg';,end
   
   if ~(strcmp(varargin{i},'p')|strcmp(varargin{i},'t')|strcmp(varargin{i},'flipy')|strcmp(varargin{i},'spm')| ...
        strcmp(varargin{i},'cor')|strcmp(varargin{i},'sagg')|strcmp(varargin{i},'ax')) 
      if isstr(varargin{i})
         basename = varargin{i}; disp(['Loading image ' basename])
      elseif ndims(varargin{i}) == 3
         array = varargin{i}; basename = 'null';
      end
   end
end

return


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ %
function [rows,cols,skip,whichslices] = getplotdims(array,print)

	if size(array,3) == 28
		rows = 6; cols = 5;skip = 1;whichslices = 1:size(array,3);
	elseif size(array,3) == 35
		rows = 6; cols = 6;skip = 1;whichslices = 1:size(array,3);
	else
		if size(array,3) > 100
			skip = floor(size(array,3) / 28);
			if print,disp(['selecting only every ' num2str(skip) 'th/rd slice.']),end
      		whichslices = 1:skip:size(array,3);
			rows = 1;cols = 1;index = 1;
			while rows * cols < ceil(size(array,3)/skip) + 1
				if mod(index,2) == 0,cols = cols + 1;
				else rows = rows + 1;
				end
				index = index + 1;
			end			
		else 
			skip = 1;
			rows = 1;cols = 1;index = 1;
			while rows * cols < ceil(size(array,3)/skip) + 1
				if mod(index,2) == 0,cols = cols + 1;
				else rows = rows + 1;
				end
				index = index + 1;
			end
			whichslices = 1:size(array,3);index = 1;
		end
	end	
return


	
