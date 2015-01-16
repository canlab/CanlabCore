function [ts,vols,chunksize] = timeseries4(coords,P,varargin)
% function [ts,vols,chunksize] = timeseries4(coords,P,[chunksize],[nochk])
%
% Simple extraction from images named in str mtx P or vols
% from voxel coordinates (not mm!) listed in coords
%
% P can be filenames in str matrix (char array)
% or 4-D array of all volume info (vols)
% (i.e., put vols in output back in as P)
%
% ts is timeseries, with fields avg and indiv for average cluster
% and individual voxels
%
% Loads images 'chunksize' at a time; default is based on memory
% size of 2^29 bytes; empty uses default
% Optional 4th argument suppresses data validity checking
%
% vols is 4-D array of all data, [x y z time(image)]
%
% Uses spm_vol and spm_read_vols and spm_slice_vol
%
% Tor Wager, 2/9/05     change from timeseries3: uses slice-by-slice
%                       method, faster.
% Tor Wager, 12/11/05   speedup for getdata with large n. voxels; cosmetic
%                       changes to output of volume method.

% -------------------------------------------------------------------
% * set up input arguments
% -------------------------------------------------------------------
maxmem = 2^28;      % note:tested the G5s up to 200 imgs, no slowdown, so no need to chunk...
chunksize = [];
    
global defaults
if isempty(defaults), spm_defaults, end

if ischar(P),
    fprintf(1,'Map vols: '); t1= clock;
    V = spm_vol(P); 
    nimages = length(V);   % number of images in data
    fprintf(1,'%3.0f s.\n', etime(clock,t1));
        
elseif ismatrix(P), 
    chunksize = NaN;
    nimages = size(P,4);    % number of images in data
else
    error('P input must be string or data matrix.');
end

if length(varargin) > 0, 
    chunksize = varargin{1};
end


if size(coords,2) ~= 3, coords = coords';, end


% -------------------------------------------------------------------
% * get extraction method
% -------------------------------------------------------------------
whslices = unique(coords(:,3));     % which slices to extract
nslices = length(whslices);         % number of z slices to extract from

extype = 'slice';
if ~isstr(P),extype = 'volume';,end
if nslices > 15, extype = 'volume';,end

% -------------------------------------------------------------------
% * read images
% -------------------------------------------------------------------

switch extype

% * slice loading method
% -------------------------------------------------------------------
case 'slice'
    
    ts.indiv = NaN .* zeros(nimages,size(coords,1)); % placeholder for data extracted
    
    for i = 1:nslices
        
        sliceno = whslices(i);               % slice number
        
    
        whcoords = find(coords(:,3) == sliceno);    % indices of in-slice coordinates 
        slcoords = coords(whcoords,:);              % coordinates to extract data from for this slice
        
        sl = timeseries_extract_slice(V,sliceno);
        
        for c = 1:size(slcoords,1)

            ts.indiv(:,whcoords(c)) = sl(slcoords(c,1),slcoords(c,2),:);
            
        end
        
        ts.avg = nanmean(ts.indiv')';
        vols = sl;       
        
end

    
% * whole-brain loading method
% -------------------------------------------------------------------
case 'volume'
        
    
if isempty(chunksize) 
    
    v = spm_read_vols(V(1));
    
    tmp = whos('v'); imgsize = tmp.bytes;
    chunksize =   floor(maxmem ./ imgsize);  %round(maxmem ./ (size(coords,1)^3 .* 16));   % images to load at a time

end

if ischar(P)

    if length(V) < chunksize
        fprintf(1,'\tChunking %3.0f into %3.0f imgs: ',length(V),chunksize);
        t1 = clock;
        fprintf(1,'Load: ');
        vols = spm_read_vols(V);
        fprintf(1,'%3.0f s. Cluster. ', etime(clock,t1));
        t1 = clock;
        [ts.indiv,ts.avg] = getdata(vols,coords,maxmem);
        fprintf(1,'%3.0f s. ', etime(clock,t1));
        
    else
        % chunk it!  load chunksize images as a whole, extract, and
        % concatenate
        ts.indiv = []; ts.avg = [];
        ind = 1;
        for i = 1:chunksize:length(V)
            t1 = clock;
            fprintf(1,'Load: %3.0f ',i);
            e = min(i+chunksize-1,length(V));
            wh = i:e;
            vols = spm_read_vols(V(wh));
            fprintf(1,'%3.0f s. Cluster. ', etime(clock,t1));
            t1 = clock;
            [indiv{ind},avg{ind}] = getdata(vols,coords,maxmem);
            ind = ind + 1;
            fprintf(1,'%3.0f s. ', etime(clock,t1));
        end
        fprintf(1,'Cat. ')
        ts.indiv = cat(1,indiv{:});
        ts.avg = cat(1,avg{:});
        clear indiv; clear avg;
    end
    
else
    % volumes already loaded
    vols = P; P = 1;
    [ts.indiv,ts.avg] = getdata(vols,coords,maxmem);
end


end         %   end switch extraction type


% -------------------------------------------------------------------
% * check for proper extraction against spm_read_vols
% -------------------------------------------------------------------

if length(varargin) > 1 | (~ischar(P)),
    % already loaded or suppress checking.
else
    fprintf(1,'Chk.\n')
    chk = check_timeseries_vals(V,ts.indiv,coords);
    if chk, keyboard, end
end

return










function [ind,avg] = getdata(vols,coords,maxmem)

    if size(coords,1) == 1  % only one voxel
        co = 1;
        ind(:,co) = squeeze(vols(coords(co,1),coords(co,2),coords(co,3),:));
        
        
    elseif size(coords,1)^3 < inf  % always do this.
        
        % time increases linearly with size of matrix; so do it in chunks.
        csz = round(sqrt(size(coords,1)));  % optimal chunk size to keep arrays as small as possible.
        indx = 1;
        for i = 1:csz:size(coords,1)
            tmp = [];
            for co = i:min(i+csz-1,size(coords,1))
                %t1 = clock;
                tmp = [tmp squeeze(vols(coords(co,1),coords(co,2),coords(co,3),:))];
                %et(co) = etime(clock,t1);
            end
            ind{indx} = tmp;
            indx = indx + 1;
            %ind = [ind squeeze(vols(coords(co,1),coords(co,2),coords(co,3),:))];
            
        end
        ind = cat(2,ind{:});
        
        
        
    else    % not a good idea speed-wise, apparently.
        tmp =  vols(coords(:,1),coords(:,2),coords(:,3),:); % the values of interest are on the 3-D diagonals of tmp
        s = size(coords,1);      % we can get the index values for diagonals by skipping elements of sv, below
        i = 1:s.^2 + s + 1:s^3;  % same as t1 = [1:size(coords,1)]'; i = sub2ind(size(tmp),t1,t1,t1)
    
        sv = (1:s^3:prod(size(tmp))) - 1; % starting values for each volume (minus one, so we add this to i)
        sv = repmat(sv,s,1) + repmat(i',1,size(sv,2));  % get the matrix of index values for each voxel at each time
        ind = tmp(sv)';
    end
    
    if size(ind,2) > 1,
        avg = nanmean(ind')';
    else
        avg = ind;
    end
    
return
    
    


        function sl = timeseries_extract_slice(V,sliceno);
        % function sl = timeseries_extract_slice(V,sliceno)
        %
        % For a given set of image names or memory mapped volumes (V)
        % extracts data from slice # sliceno and returns an X x Y x time
        % matrix of data.
        % uses spm_slice_vol.m
      
            if isstr(V), V = spm_vol(V);,end
            
            mat = spm_matrix([0 0 sliceno]);     % matrix for spm_slice_vol
        
            for i = 1:length(V)
                sl(:,:,i) = spm_slice_vol(V(i),mat,V(i).dim(1:2),0); 
            end
        
        return
        
        
        

function chk = check_timeseries_vals(V,dat,coords)
% function chk = check_timeseries_vals(V,dat,coords)
%
% Checks a random subset of up to 5 images against extracted data values
% to make sure the right data is being extracted.
%
% V is memory-mapped volumes (see spm_vol)
% dat is extracted data
% coords is voxel coordinates, n rows x 3 columns
%
% tor wager

n = min(5,length(V));

% get random set of n images, up to 5
wh = randperm(length(V));
wh = wh(1:n);

% select these rows in V and dat
% dat is [images, coordinates]
V = V(wh);
dat = dat(wh,:);

% get random set of nn coordinates, up to 5

nc = size(coords,1);
nn = min(5,nc);
whc = randperm(nc);
whc = whc(1:nn);
coords = coords(whc,:);

% select these columns in dat
dat = dat(:,whc);

v = spm_read_vols(V);
for i = 1:nn    % for each coordinate
    dat2(:,i) = squeeze(v(coords(i,1),coords(i,2),coords(i,3),:));
end

chk = dat - dat2;
chk = any(chk(:));

if chk, warning('Problem with timeseries3!! Extracted data do not match expected values.  Quitting at error');,end

return


    
