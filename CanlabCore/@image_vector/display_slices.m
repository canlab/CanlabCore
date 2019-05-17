function whsl = display_slices(dat, varargin)
% Creates 3 separate montage views - ax, cor, sagg in a special figure window
%
% - By default, a figure with axial, coronal, and saggital montages are
% created.  But it's also easy to create nicer-looking separate figures with
% only one view.
% - Also easy to select slices: Enter 'startslice' and 'endslice', each followed by values, to display
% a range of slices. Default is entry in mm coordinate, but 'voxelspace' allows you to specify them in 
% voxel coodinates, i.e., in slices in the actual data image.
%
% - This function is good for getting a quick view of the contents of an
% image, e.g., for structural images.
% - The montage.m method and fmridisplay are good for displaying blobs on
% standard brains, and are much slower.
%
% :Usage:
% ::
%
%    display_slices(dat, [myview], ['spacing', slicespacing], ['vertical'])
%
% :Examples:
% ::
%
%     display_slices(dat);
%     display_slices(dat, 'coronal');
%     display_slices(dat, 'saggital', 'spacing', 10); % 10 mm spacing
%     display_slices(dat, 'saggital', 'spacing', 10, 'vertical');
%     display_slices(dat, 'saggital', 'slices_per_row', 12);
%     display_slices(dat, 'axial', 'slices_per_row', 20, 'spacing', 4, 'startslice', -10, 'endslice', 10); % z = -10 to 10 mm, 4 mm spacing
%
% Make a montage of a subset of slices and make it fill the figure:
% create_figure('sagg')
% display_slices(dat, 'saggital', 'slices_per_row', 4, 'spacing', 3, 'startslice', -20, 'endslice', 20)
% set(gca, 'Position', [.05 .05 .9 1]);
% colormap gray
%
% ..
%    Tor Wager, Aug 2012. Update: Jan 2018
% ..

do3views = 1;     % composite with all 3 views (default) in special figure
myview = 'axial'; % default for single view
spacing = 1;
dovertical = false;
stackorient = 2;  % 2 for horiz, 1 for vertical
dotight = true;   % tight: show only non-zero slices
s = 12;           % slices per row
% startslice = [];  % first slice to show; set later 
% endslice = [];    % last slice to show (empty = all); set later 
entered_mm_coords = true;  % interpret entries in mm or voxel space coordinates

% Process inputs - single case
% in current axes
% ------------------------------------------------------------------------
for i = 1:length(varargin)
    do3views = false;  % if any optional arguments entered, turn off 3 views
    
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case {'axial', 'saggital', 'sagittal', 'coronal'}, myview = varargin{i};
            case 'spacing', spacing = varargin{i + 1};
            case 'vertical', stackorient = 1;
            case 'slices_per_row', s = varargin{i + 1};
            case 'startslice', startslice = varargin{i + 1}; 
            case 'endslice', endslice = varargin{i + 1}; 
            case 'mm', entered_mm_coords = true;            % we have entered coords in mm [default]
            case 'voxelspace', entered_mm_coords = false;   % entered values in voxel space
            otherwise
                disp('Warning: Unknown input.')
        end
    end
end

% wh = strcmp(varargin, 'voxelspace'); % entered in voxel space
% if any(wh), entered_mm_coords = false; end


% ------------------------------------------------------------------------
% Multi-view mode (default)
%
if do3views
    fh = create_figure('slice_display');
    axis off
    ax1 = axes('OuterPosition', [.05 .05 .20 .90], 'Position', [.05 .05 .20 .90]);
    display_slices(dat, 'saggital', 'spacing', 8, 'vertical');
    
    ax2 = axes('OuterPosition', [.229 .05 .8 .150]); %[.20 .05 .8 .30]);
    set(ax2, 'Position', get(ax2, 'OuterPosition'));
    display_slices(dat, 'coronal', 'spacing', 8);
    
    ax3 = axes('OuterPosition', [.20 .25 .8 .70]);
    set(ax3, 'Position', get(ax3, 'OuterPosition'));
    display_slices(dat, 'axial');
    
    colormap gray
    return
end

% ------------------------------------------------------------------------
% Main function for single montage
%

if stackorient == 2, platestack = 1; else platestack = 1; end

vdat = reconstruct_image(dat);

% flip x dim if x voxel size is negative. other display/write methods also
% do this.
if sign(dat.volInfo.mat(1)) < 0
    vdat = flipdim(vdat, 1);
end

% Find out which dimension to vary slices along
% ------------------------------------------------------------------------
switch myview
    case 'axial'
        wh_col = 3;
    case {'saggital', 'sagittal'}
        wh_col = 1;
    case 'coronal'
        wh_col = 2;
    otherwise error('unknown slice view.')
end

% Get vector of slices in mm or voxel space, depending on input args
% ------------------------------------------------------------------------
% re-set default slices if we have not entered them
% depends on voxel or mm. if mm, convert to voxels.
if ~any(strcmp(varargin, 'startslice'))
    if entered_mm_coords  % Figure out mm coords of first slice
        xyzvoxel = voxel2mm([1 1 1], dat.volInfo.mat);
        startslice = xyzvoxel(wh_col); 
    else
        startslice = 1;
    end
end

if ~any(strcmp(varargin, 'endslice'))
    if entered_mm_coords % Figure out mm coords of last slice
        xyzvoxel = voxel2mm(size(vdat)', dat.volInfo.mat);
        endslice = xyzvoxel(wh_col);
    else
        endslice = size(vdat, wh_col);
    end
end

whsl = startslice:spacing:endslice;  % voxel or mm


% mm to voxel conversion
% ------------------------------------------------------------------------
if entered_mm_coords
    % Convert from mm 2 voxel, if we have entered mm coordinates
    
    xyzmm = zeros(length(whsl), 3);
    
    xyzmm(:, wh_col) = whsl;
    
    xyzvoxel = mm2voxel(xyzmm', dat.volInfo.mat, 1);
    whsl = xyzvoxel(:, wh_col)';
    
end

% Now everything is in voxel coords from here on out.

if dotight
   % figure out which slices to skip because they are empty.
   switch myview
    case 'axial'
        subvol = abs(vdat(:, :, whsl));
        sumabsval = sum(sum(subvol, 1), 2);
        
    case {'saggital', 'sagittal'}
        subvol = abs(vdat(whsl, :, :));
        sumabsval = sum(sum(subvol, 2), 3);
        
    case 'coronal'
        subvol = abs(vdat(:, whsl, :));
        sumabsval = sum(sum(subvol, 1), 3);
        
    otherwise error('unknown slice view.')
   end

   wh_omit = sumabsval < 100 * eps;
   whsl(wh_omit) = []; 

end

% Eliminate empty areas of image (skip dim we're varying along)
% ------------------------------------------------------------------------
if dotight
    vdat = eliminate_empty_areas(vdat, myview);
end

% stack images to display together
% ------------------------------------------------------------------------

% starting and ending values for slices
% indices into whsl
nrows = ceil(length(whsl) ./ s);

for j = 1:nrows
    rowslices{j} = whsl((j-1)*s+1 : min(j*s, length(whsl)));
end

plate = cell(nrows, 1);

for j = 1:nrows
    stacked = [];
    
    for i = rowslices{j}  %st(j):en(j)
        
        switch myview
            case 'axial'
                stacked{i} = vdat(:, :, i);
                stacked{i} = rot90(stacked{i});
            case {'saggital', 'sagittal'}
                stacked{i} = squeeze(vdat(i, :, :));
                stacked{i} = rot90(stacked{i});
            case 'coronal'
                stacked{i} = squeeze(vdat(:, i, :));
                stacked{i} = rot90(stacked{i});
            otherwise error('unknown slice view.')
        end
        
    end
    
    stacked = cat(stackorient, stacked{:});
    
    plate{j} = stacked;
    
end % rows/plates

maxrows = max(cellfun(@(x) size(x, 1), plate, 'UniformOutput', true));
maxcols = max(cellfun(@(x) size(x, 2), plate, 'UniformOutput', true));
padval = 0; % mean(plate{j}(:)); % value to use for pad

% Pad if needed
for j = 1:nrows
    if size(plate{j}, 2) < maxcols
        % pad
        plate{j} = [plate{j} padval .* ones(maxrows, maxcols - size(plate{j}, 2))];
    end
    
end % rows

plate = cat(platestack, plate{:});

% create image
% ------------------------------------------------------------------------
imagesc(plate)
axis image
set(gca, 'YDir', 'Reverse')
axis off

end

% function vdat = eliminate_empty_areas(vdat)
% % Eliminate empty areas of image
% % ------------------------------------------------------------------------
% % this will mess up conversion to voxel coords from mm if empty slices!
% % turn off in direction of variation
% 
%     nullvox = vdat == 0 | isnan(vdat);
%     bottom = all(nullvox, 3);
%     
%         nullx = squeeze(all(bottom, 2));
%         vdat(nullx, :, :) = [];
%     
%         nully = squeeze(all(bottom, 1));
%         vdat(:, nully, :) = [];
%     
%         side = squeeze(all(nullvox, 1));
%         nullz = squeeze(all(side, 1));
%         vdat(:, :, nullz) = [];
%         
% end

function vdat = eliminate_empty_areas(vdat, myview)
% Eliminate empty areas of image
% ------------------------------------------------------------------------
% this will mess up conversion to voxel coords from mm if empty slices!
% turn off in direction of variation

    nullvox = vdat == 0 | isnan(vdat);
    bottom = all(nullvox, 3);
    
    if strcmp(myview, 'saggital') | strcmp(myview, 'sagittal')
        % skip
    else
        nullx = squeeze(all(bottom, 2));
        vdat(nullx, :, :) = [];
    end
    
    if strcmp(myview, 'coronal')
        % skip
    else
        nully = squeeze(all(bottom, 1));
        vdat(:, nully, :) = [];
    end
    
    if strcmp(myview, 'axial')
        % skip
    else
        side = squeeze(all(nullvox, 1));
        nullz = squeeze(all(side, 1));
        vdat(:, :, nullz) = [];
    end
    
end
