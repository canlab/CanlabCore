function cl = orthviews(image_obj, varargin)
% Orthviews display (SPM) for CANlab image_vector (or fmri_data, statistic_image) object
%
% Usage:
% orthviews(image_obj, varargin)
%
% Options
%   'posneg' input generates orthviews using solid colors.
%   'largest_region' to center the orthviews on the largest region in the
%   image
%
% Copyright Tor Wager, 2011
%
% ---------------------------------------------------------------

overlay = which('SPM8_colin27T1_seg.img');

doposneg = 0;
doreg = 0;
input_handle=[];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'posneg',  doposneg = 1;
                
                % functional commands
            case 'overlay', overlay = varargin{i + 1}; varargin{i + 1} = [];
                
            case {'han', 'handle', 'input_handle'}, input_handle = varargin{i+1};
                
            case {'largest_region', 'largest_cluster'}, doreg = 1;
                
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


for i = handle_indices
    
    cl{i} = iimg_indx2clusters(image_obj.dat(:, min(i, size(image_obj.dat,2))), image_obj.volInfo); %min() is to differentiate between which subplot to plot on and which image to plot.  used when this is called with orthviews_multiple_objs.  Yoni 11/14
    
    if doposneg
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
    
    spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 1 1], [.5 .5 .5], [1 .5 0]);
    
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



