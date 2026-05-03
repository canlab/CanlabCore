function cl = orthviews(image_obj, varargin)
% orthviews Orthviews display (SPM) for a statistic_image object.
%
% Displays each image in a statistic_image as an SPM orthviews panel,
% honoring the .sig field so only suprathreshold voxels are shown.
%
% :Usage:
% ::
%
%     cl = orthviews(image_object)
%     cl = orthviews(image_object, handle_number of existing orthviews)
%
% :Features:
%
%   - Uses cluster_orthviews.m
%   - Will take inputs to cluster_orthviews
%   - e.g., orthviews(pstat, 'unique', 'solid');
%
% Output is clusters structure (see also region.m).
%
% :Inputs:
%
%   **image_obj:**
%        A statistic_image object. Each column of .dat is rendered in a
%        separate orthviews panel; .sig is used as a mask if non-empty.
%
% :Optional Inputs:
%
%   **'handle' / 'han' / 'input_handle':**
%        Followed by the handle number(s) of an existing orthviews
%        figure to display into instead of opening a new one.
%
%   **'largest_region':**
%        Center the orthviews crosshair on the largest region in the
%        first image.
%
%   **'overlay':**
%        Followed by a filename to use as the underlay image. Default:
%        which('fmriprep20_template.nii.gz').
%
%   Additional optional arguments are passed through to cluster_orthviews.
%
% :Outputs:
%
%   **cl:**
%        Cell array of clusters structures, one per displayed image
%        (see also region.m).
%
% :Examples:
% ::
%
%    % T-test, Construct a stats_image object, threshold and display:
%    statsimg = ttest(fmridat, .001, 'unc');
%
%    % Re-threshold and display:
%    statsimg = threshold(statsimg, .000001, 'unc');
%    orthviews(statsimg);
%
%    statsimg = threshold(statsimg, .01, 'fdr');
%    orthviews(statsimg);
%
%    % Create an orthviews and view at multiple thresholds in different panes:
%    overlay = which('SPM8_colin27T1_seg.img');
%    spm_check_registration(repmat(overlay, n, 1));
%    statsimg = ttest(fmridat);
%    statsimg = threshold(statsimg, .001, 'unc');
%    orthviews(statsimg, 'handle', 1);
%
%    statsimg = threshold(statsimg, .000001, 'unc');
%    orthviews(statsimg, 'handle', 2);
%
% :See also:
%   - statistic_image.multi_threshold
%   - cluster_orthviews
%   - region

input_handle = [];
cl = [];
doreg=0;
overlay = which('fmriprep20_template.nii.gz');
%dounique = 0; uniquestr = 'nounique';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'han', 'handle', 'input_handle'}, input_handle = varargin{i+1};
            case 'largest_region', doreg = 1;  
            case 'overlay', overlay = varargin{i + 1}; varargin{i + 1} = [];    
            %case 'unique', dounique = 1; uniquestr = 'unique';
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(input_handle)
    % re-initialize new orthviews
    
    n = size(image_obj.dat, 2);

    spm_check_registration(repmat(overlay, n, 1));
    
    handle_indices = 1:n;
else
    
    % use existing
    handle_indices = input_handle;
    
end

image_indx = 1;

if isempty(image_obj.sig), image_obj.sig = true(size(image_obj.dat)); end

% replace missing voxels if necessary
image_obj.sig = zeroinsert(image_obj.removed_voxels, image_obj.sig);
image_obj.dat = zeroinsert(image_obj.removed_voxels, image_obj.dat);

% do not do this, because it may replace deliberately removed imgs
% image_obj = replace_empty(image_obj);


for i = handle_indices
    
    if size(image_obj.sig, 2) < image_indx
        disp('No sig vector; displaying all voxels')
    end
    
    % use .sig field if we have it, otherwise use dat
    if size(image_obj.sig, 2) >= image_indx && any(image_obj.sig(:, image_indx))  
        cl{i} = iimg_indx2clusters(image_obj.dat(:, image_indx) .* double(image_obj.sig(:, image_indx)), image_obj.volInfo);
    else
        cl{i} = iimg_indx2clusters(image_obj.dat(:, image_indx), image_obj.volInfo);
        %disp('No sig vector; displaying all voxels')
    end
    
    if ~isempty(cl{i})
        if ~isempty(image_obj.image_names) && ischar(image_obj.image_names)
            
            if size(image_obj.image_names, 1) < image_indx
                image_obj.image_names = char(image_obj.image_names, 'NO NAME');
            end
            
            fprintf('Displaying image %3.0f, %3.0f voxels: %s\n', i, sum(cat(1, cl{i}.numVox)), image_obj.image_names(image_indx, :));
        end
        
        cluster_orthviews(cl{i}, 'add', 'handle', i, varargin{:});
                
    elseif size(image_obj.image_names, 1) >= image_indx
        fprintf('Image %3.0f empty: %s\n', i, image_obj.image_names(image_indx, :));
    else 
        fprintf('Image %3.0f empty\n', i);
    end
    
    image_indx = image_indx + 1;
end

%if ~dounique
    spm_orthviews_change_colormap([0 0 1], [1 1 0], [0 1 1], [.5 .5 .5], [1 .5 0]);
%end

% set crosshairs
if doreg
    [~,wh] = max(cat(1,cl{1}.numVox));
    spm_orthviews('Reposition', cl{1}(wh).mm_center);
else
    spm_orthviews('Reposition', [0 0 0]); 
end
drawnow

% Set colormap
spm_orthviews_hotcool_colormap(image_obj.dat(:), 0);

end % function


