function o = slices(obj, varargin)
% Create a montage of single-slice results for every image in an image_vector object
%
% :Usage:
% ::
%
%    o = slices(obj, 'orientation', [orientation], 'slice', [slice_mm], 'nimages', [nimgs])
%
% obj is an image_vector, fmri_data, or statistic_image object with
% multiple images (only the first 64 will display), which are stored as
% columns in its .dat field.
%
% :Optional Inputs:
%
%   **orientation:**
%        can be followed by 'saggital', 'axial', or 'coronal'
%
%   **slice_mm:**
%        is followed by the mm coord of the slice to display; default = 0
%
%   **nimgs:**
%        can be followed by the number of images to display, 1:nimgs
%
%   **names:**
%        is followed by a cell array of names for the images.
%
%   **color:**
%        is followed by color vector or string specification. default is
%        color-mapped with split colors (hot/cool) for pos and neg effects.
%
%   **outline:**
%        is followed by a color vector for outline around blobs.
%
% The output, o, is an fmridisplay object.
%
% This function uses fmridisplay objects, and may be memory-intensive for
% older computers.
%
% *Common Errors:*
%
% This function uses the volInfo.cluster field. If you create a mask in an
% ad hoc way, this field may not be updated.  use this to fix:
%   - mask = reparse_contiguous(mask);
%
% :Examples:
% ::
%
%    slices(dat);
%    slices(dat, 'orientation', 'axial');
%    slices(dat, 'slice', -5);                 % display sagg at x = -5
%    o = slices(dat, 'names', terms); % use 'terms' var as names
%
%    o2 = slices(all_chi2_images, 'orientation', 'saggital', 'slice', 0);
%
% ..
%    Copyright 2011, Tor Wager
% ..

slice_mm = 0;
my_orientation = 'saggital';
nimgs = size(obj.dat, 2);
dosplitcolor = 1;
outlinecolor = [];
for i = 1:nimgs, names{i, 1} = sprintf('Img %3.0f', i); end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'orientation', my_orientation = varargin{i+1}; varargin{i+1} = [];
            case 'slice', slice_mm = varargin{i+1};
            case 'nimages', nimgs = varargin{i+1};
            case 'names', names = varargin{i+1};
                
            case 'color', dosplitcolor = 0; color = varargin{i + 1};
            case 'outline', outlinecolor = varargin{i + 1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

numax = [ceil(sqrt(nimgs)) floor(sqrt(nimgs))];
if prod(numax) < nimgs, numax(2) = numax(2) + 1; end

if nimgs > 64, disp('Displaying only first 64 images!'); end

m = mean(obj);

o = fmridisplay;
f1 = create_figure('slice_montage', numax(1), numax(2));
colormap gray

for i = 1:prod(numax)
    newaxhan(i) = subplot(numax(1), numax(2), i);
    axis off;
end



% For each slice
for i = 1:nimgs
    
    %o = removeblobs(o);
    o = montage(o, my_orientation, 'slice_range', [slice_mm slice_mm]);
    enlarge_axes(gcf, .8)
    
    % add blobs
    m.dat = obj.dat(:, i);
    
    if isa(obj, 'statistic_image')
        m.dat(~obj.sig(:, i)) = 0;
    end
    
    m = reparse_contiguous(m);
    
    cl = region(m);
    
    %     cluster_orthviews(cl, {[1 0 0]}, 'solid');
    %     %spm_orthviews_name_axis(names{i}, i);
    %     spm_orthviews('Position', [0 0 0]);
    %     s = input('Press a key');
    %
    if dosplitcolor
        o = addblobs(o, cl, 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}, 'wh_montages', i);
    else
        o = addblobs(o, cl, 'color', color, 'wh_montages', i);
    end
    
    if ~isempty(outlinecolor)
        o = addblobs(o, cl, 'color', outlinecolor, 'wh_montages', i, 'outline');
    end
    
    title(names{i}, 'FontSize', 24)
    
    % copy to composite figure
    newhan(i) = copyobj(o.montage{i}.axis_handles(1), f1);
    
    o.montage{i}.axis_handles = newhan(i);
    
    % close the original
    close
    
    newpos = get(newaxhan(i), 'Position');
    set(newhan(i), 'Position', newpos);
    drawnow
    
end


end % function
