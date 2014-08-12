function fastmontage(dat, varargin)
% fastmontage(dat, [myview], ['spacing', slicespacing], ['vertical'])
%
% fastmontage(dat);  Creates 3 separate montage views - ax, cor, sagg
%                    In special figure window
%
% fastmontage(dat, 'coronal');
% fastmontage(dat, 'saggital', 'spacing', 10);
% fastmontage(dat, 'saggital', 'spacing', 10, 'vertical');
% fastmontage(dat, 'saggital', 'slices_per_row', 12);
%
% Tor Wager, Aug 2012

do3views = 1;     % composite with all 3 views (default) in special figure
myview = 'axial'; % default for single view
spacing = 1;
dovertical = 0;
stackorient = 2;  % 2 for horiz, 1 for vertical
dotight = 1;      % tight: show only non-zero slices
s = 12;           % slices per row

% Process inputs - single case
% in current axes
% ------------------------------------------------------------------------
for i = 1:length(varargin)
    do3views = 0;  % if any optional arguments entered, turn off 3 views
    
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case {'axial', 'saggital', 'coronal'}, myview = varargin{i};
            case 'spacing', spacing = varargin{i + 1};
            case 'vertical', stackorient = 1;
            case 'slices_per_row', s = varargin{i + 1};
                
            otherwise
                disp('Warning: Unknown input.')
        end
    end
end


% ------------------------------------------------------------------------
% Multi-view mode (default)
%
if do3views
    fh = create_figure('fastmontage');
    ax1 = axes('OuterPosition', [.05 .05 .20 .90], 'Position', [.05 .05 .20 .90]);
    fastmontage(dat, 'saggital', 'spacing', 8, 'vertical');
    
    ax2 = axes('OuterPosition', [.229 .05 .8 .150]); %[.20 .05 .8 .30]);
    set(ax2, 'Position', get(ax2, 'OuterPosition'));
    fastmontage(dat, 'coronal', 'spacing', 8);
    
    ax3 = axes('OuterPosition', [.20 .25 .8 .70]);
    set(ax3, 'Position', get(ax3, 'OuterPosition'));
    fastmontage(dat, 'axial');
    
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

if dotight
    nullvox = vdat == 0 | isnan(vdat);
    bottom = all(nullvox, 3);
    nullx = squeeze(all(bottom, 2));
    vdat(nullx, :, :) = [];
    
    nully = squeeze(all(bottom, 1));
    vdat(:, nully, :) = [];
    
    side = squeeze(all(nullvox, 1));
    nullz = squeeze(all(side, 1));
    vdat(:, :, nullz) = [];
    
end

% Spacing and which slices
switch myview
    case 'axial'
        z = size(vdat, 3);
    case 'saggital'
        z = size(vdat, 1);
    case 'coronal'
        z = size(vdat, 2);
    otherwise error('unknown slice view.')
end
whsl = 1:spacing:z;

nrows = ceil(length(whsl) ./ s);

% starting and ending values for slices
%stepby = ceil(length(whsl)/nrows);
% indices into whsl
for j = 1:nrows
    rowslices{j} = whsl((j-1)*s+1 : min(j*s, length(whsl))); 
end
% st = whsl(1:8:length(whsl));
% en = min(st+stepby-1, length(whsl));

plate = cell(nrows, 1);

for j = 1:nrows
    stacked = [];
    
    for i = rowslices{j}  %st(j):en(j)
        
        switch myview
            case 'axial'
                stacked{i} = vdat(:, :, i);
                stacked{i} = rot90(stacked{i});
            case 'saggital'
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
    
    if j > 1 && size(plate{j}, 2) < size(plate{1}, 2)
        % pad
        fullsize = size(plate{1});
        padval = 0; % mean(plate{j}(:)); % value to use for pad
        plate{j} = [plate{j} padval .* ones(fullsize(1), fullsize(2) - size(plate{j}, 2))];
    end
    
end % rows


plate = cat(platestack, plate{:});

imagesc(plate)
axis image
set(gca, 'YDir', 'Reverse')
axis off

end
