function orthviews(image_obj, varargin)
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

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'posneg',  doposneg = 1;
                
                % functional commands
            case 'overlay', overlay = varargin{i + 1}; varargin{i + 1} = [];
                
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

% replace missing voxels if necessary
image_obj = replace_empty(image_obj);

n = size(image_obj.dat, 2); % images are columns here

if n > 24
    disp('Warning! Only first 24 images can displayed with spm_orthviews.');
    n = 24;
end


spm_check_registration(repmat(overlay, n, 1));

for i = 1:n
    
    cl{i} = iimg_indx2clusters(image_obj.dat(:, i), image_obj.volInfo);
    
    if doposneg
        warning off
        cl{i} = cluster2region(cl{i});
        warning on % name field warning
        [clpos{i}, clneg{i}] = posneg_separate(cl{i}, 'Z');
        
        if ~isempty(clpos{i})
            cluster_orthviews(clpos{i}, {[1 1 0]}, 'add', 'handle', i, 'solid');
        end
        
        if ~isempty(clneg{i})
            cluster_orthviews(clneg{i}, {[0 0 1]}, 'add', 'handle', i, 'solid');
        end
        
    else
        cluster_orthviews(cl{i}, 'add', 'handle', i);
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



