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

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'posneg',  doposneg = 1;
                
                % functional commands
            case 'overlay', overlay = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'han', 'handle', 'input_handle'}, input_handle = varargin{i+1};
                
            case {'largest_region', 'largest_cluster'}, doreg = 1;
                
            case 'unique'
                dounique = true;
                
            case 'continuous'
                force_continuous = true;
                
            case 'parcels'
                doparcels = true; dounique = true;
                
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

% autodetect parcels and turn on 'unique' for maps with few values 
% (default = do this)
if ~doparcels && ~force_continuous
    nvalues = length(unique(image_obj.dat(:)));
    
    if nvalues < 300, dounique = true; doparcels = true; end
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
            cluster_orthviews(cl{i}, 'add', 'handle', i);
        else
            fprintf('Image %3.0f empty\n', i);
        end
    end
    
    %spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 1 1], [.5 .5 .5], [1 .5 0]);
    spm_orthviews_change_colormap([.5 0 1], [1 1 0]);
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



