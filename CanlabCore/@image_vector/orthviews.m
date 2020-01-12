function cl = orthviews(image_obj, varargin)
% Orthviews display (SPM) for CANlab image_vector (or fmri_data, statistic_image) object
%
% *Usage:
% ::
%
%    orthviews(image_obj, varargin)
%
% :Optional Inputs:
%   **'posneg':**
%        input generates orthviews using solid colors, separated for positive- and negative-valued voxels.
%
%   **'largest_region':**
%        to center the orthviews on the largest region in the image
%
%   **'overlay':**
%        followed by name of anatomical image to use as underlay
%
%   **'unique':**
%        plot groups of contiguous voxels in different, unique colors
%        - autodetect feature turns this on when there are <300 unique
%        values, or for 'atlas' objects. See 'continuous' below to force continuous.
%
%   **'parcels':**
%        plot voxels with each unique value in a different color
%        useful for images in which values indicate parcel or network
%        membership
%        Default is 'unique' for <300 values, 'continuous' for more.  Use
%        'unique' to force unique-valued colors.
%
%   **'continuous':**
%        Plot voxels color-mapped according to their values.
%        Default is 'unique' for <300 values, 'continuous' for more.  Use
%        'unique' to force unique-valued colors.
%
%   **'input_handle':**
%        followd by integer for which orthviews display to add to.
%        spm_orthviews can display multiple linked orthviews sections. If
%        you have initialized a display with multiple images, Enter 1 to
%        add to the first one, 2 for the second, and so forth.
%      
%   **'add':**
%        Add to the first orthviews display (e.g., 'input_handle', 1)
%
%   **'color':**
%       Followed by cell with 3-element rgb vector, 
%       e.g., 'color', {[1 0 1]}
%       This is superseded by 'unique' and 'posneg' options.
%
% ..
%    Copyright Tor Wager, 2011
% ..

% Aug 2017: Tor Wager edited to allow flexibility in colors; unique,
% parcel, continuous. and autodetect.

%overlay = which('SPM8_colin27T1_seg.img');
overlay = which('keuken_2014_enhanced_for_underlay.img');

doposneg = false;
doreg = false;
input_handle=[];
dounique = false;
doparcels = false;
force_continuous = false;
additional_inputs = {};

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'posneg',  doposneg = 1;
                
                % functional commands
            case 'overlay', overlay = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'han', 'handle', 'input_handle'}, input_handle = varargin{i+1};
                
            case 'add', input_handle = 1;
                
            case {'largest_region', 'largest_cluster'}, doreg = 1;
                
            case 'unique'
                dounique = true;
                
            case 'continuous'
                force_continuous = true;
                
            case 'parcels'
                doparcels = true; dounique = true;
                
            case 'color'
                additional_inputs{end + 1} = varargin{i + 1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ~exist(overlay, 'file')
    overlay = which(overlay);
end

if isempty(overlay) || ~exist(overlay, 'file')
    disp('Cannot find overlay image:')
    disp(overlay)
    return
end

if isempty(input_handle)
    % re-initialize new orthviews
    
    n = size(image_obj.dat, 2);
    if n > 24
        disp('Warning! Only first 24 images can displayed with spm_orthviews.');
        n = 24;
    end

    spm_check_registration(repmat(overlay, n, 1));
    
    handle_indices = 1:n;

else
    
    % use existing
    handle_indices = input_handle;

end

% replace missing voxels if necessary
image_obj = replace_empty(image_obj);

% enforce double, same var types across all objects
image_obj.dat = double(image_obj.dat);

% autodetect parcels and turn on 'unique' for maps with few values 
% or atlas objects.
% (default = do this)
if isa(image_obj, 'atlas') && ~force_continuous
    dounique = true; doparcels = true; 
end
    
u = unique(image_obj.dat(:));

if ~doparcels && ~force_continuous && all(u == round(u))
    % auto-select parcel unique color plot if there are few unique values and values are integers
    nvalues = length(u);
    
    if nvalues > 2 && nvalues < 300, dounique = true; doparcels = true; end
end

for i = handle_indices
    
    wh_image = min(i, size(image_obj.dat,2));  %min() is to differentiate between which subplot to plot on and which image to plot.  used when this is called with orthviews_multiple_objs.  Yoni 11/14
    
    if doparcels
        
        cl{i} = region(get_wh_image(image_obj, wh_image), 'unique_mask_values');
    
    else
        cl{i} = region(get_wh_image(image_obj, wh_image), 'contiguous_regions');
        
    end
    
    % old:
    %cl{i} = iimg_indx2clusters(image_obj.dat(:, wh_image), image_obj.volInfo); 
    
    if dounique
        
        %colors = custom_colors([1 1 0], [0 0 1], length(cl{i}));
        
        if ~isempty(cl{i})
            cluster_orthviews(cl{i}, 'unique', 'add', 'handle', i);
        else
            fprintf('Image %3.0f empty\n', i);
        end
        
    
    elseif doposneg
        
        warning off
        cl{i} = cluster2region(cl{i});
        warning on % name field warning
        
        [clpos{i}, clneg{i}] = posneg_separate(cl{i}, 'Z');
        
        if ~isempty(clpos{i})
            cluster_orthviews(clpos{i}, {[1 1 0]}, 'add', 'handle', i, 'solid');
        else
            fprintf('Image %3.0f: No positive clusters\n', i);
        end
        
        if ~isempty(clneg{i})
            cluster_orthviews(clneg{i}, {[0 0 1]}, 'add', 'handle', i, 'solid');
        else
            fprintf('Image %3.0f: No negative clusters\n', i);
        end
        
    else
        if ~isempty(cl{i})
            cluster_orthviews(cl{i}, additional_inputs{:}, 'add', 'handle', i);
        else
            fprintf('Image %3.0f empty\n', i);
        end
    end
    

end

% Set colormap
% If one image, we can use 'hot/cool' split, which shows zero point
% If multiple images, default to a single colormap, as the entire
% figure can use only one colormap.

if length(handle_indices) > 1
    %spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 1 1], [.5 .5 .5], [1 .5 0]);
    %spm_orthviews_change_colormap([.5 0 1], [1 1 0]);
    spm_orthviews_change_colormap([.2 .2 .6], [1 1 0]);  % slate to yellow
    
else % one image.
    spm_orthviews_hotcool_colormap(image_obj.dat(:), 0);
end

% center image
spm_orthviews('Reposition', [0 0 0]);

% if requested, try to center on largest region
if doreg
    [~,wh] = max(cat(1,cl{1}.numVox));
    spm_orthviews('Reposition', cl{1}(wh).mm_center);
end 
drawnow;

end % function



