function o2 = canlab_results_fmridisplay(input_activation, varargin)

% usage: function canlab_results_fmridisplay(input_activation, [optional inputs])
% Tor Wager
% 1/27/2012
% purpose:  This function display fmri results.
%
% input:    input_activation - nii, img,
%           This image has the blobs you want to
%           display. You can also enter a cl "clusters" structure or
%           "region" object.
%
%           you can also get a thresholded image like the examples used here
%           from a number of places - by thresholding your results in SPM
%           and using "write filtered" to save the image, by creating masks
%           from meta-analysis or anatomical atlases, or by using
%           mediation_brain_results, robust_results_threshold,
%           robust_results_batch_script, threshold_imgs, or object
%           oriented tools including fmri_data and statistic_image objects.
%
% optional inputs:
% -------------------------------------------------------------------------
% 'noblobs' : do not display blobs
% 'nooutline' : do not display blob outlines
% 'addmontages' : when entering existing fmridisplay obj, add new montages
% 'noremove' : do not remove current blobs when adding new ones
% 'outlinecolor : followed by new outline color
% 'splitcolor' : followed by 4-cell new split colormap colors (help fmridisplay or edit code for defaults as example)
% 'montagetype' : 'full' for full montages of axial and sagg slices. 
%                 'compact' [default] for single-figure parasagittal and
%                 axials slices.
%                 'compact2': like 'compact', but fewer axial slices.
%
% * Other inputs to addblobs (fmridisplay method) are allowed, e.g., 'cmaprange', [-2 2], 'trans'
% See help fmridisplay
% e.g., 'color', [1 0 0]
%
% You can also input an existing fmridisplay object, and it will use the
% one you have created rather than setting up the canonical slices.
%
% example script:
% -------------------------------------------------------------------------
% input_activation = 'Pick_Atlas_PAL_large.nii';
%
% % set up the anatomical underlay and display blobs
% % (see the code of this function and help fmridisplay for more examples)
%
% o2 = canlab_results_fmridisplay(input_activation);
%
% %% ========== remove those blobs and change the color ==========
%
% cl = mask2clusters(input_activation);
% removeblobs(o2);
% o2 = addblobs(o2, cl, 'color', [0 0 1]);
%
% %% ========== OR
%
% r = region(input_activation);
% o2 = removeblobs(o2);
% o2 = addblobs(o2, r, 'color', [1 0 0]);
%
% %% ========== Create empty fmridisplay object on which to add blobs:
% o2 = canlab_results_fmridisplay
%
% %% ========== If you want to start over with a new fmridisplay object,
% % make sure to clear o2, because it uses lots of memory
%
% % This image should be on your path in the "canlab_canonical_brains" subfolder:
%
% input_activation = 'pain-emotion_2s_z_val_FDR_05.img';
% clear o2
% close all
% o2 = canlab_results_fmridisplay(input_activation);
%
% %% ========== save PNGs of your images to insert into powerpoint, etc.
% % for your paper/presentation
%
% scn_export_papersetup(400);
% saveas(gcf, 'results_images/pain_meta_fmridisplay_example_sagittal.png');
%
% scn_export_papersetup(350);
% saveas(gcf, 'results_images/pain_meta_fmridisplay_example_sagittal.png');
%
% Change colors, removing old blobs and replacing with new ones:
% o2 = canlab_results_fmridisplay(d, o2, 'cmaprange', [.3 .45], 'splitcolor', {[0 0 1] [.3 0 .8] [.9 0 .5] [1 1 0]}, 'outlinecolor', [.5 0 .5]);


if ~which('fmridisplay.m')
    disp('fmridisplay is not on path.  it is in canlab tools, which must be on your path!')
    return
end

if nargin == 0
    o2 = canlab_results_fmridisplay(region(), 'noblobs', 'nooutline');
    return
end

if ischar(input_activation)
    cl = mask2clusters(input_activation);
    
elseif isstruct(input_activation) || isa(input_activation, 'region')
    cl = input_activation;
    
elseif isa(input_activation, 'image_vector')
    cl = region(input_activation);
    
else
    error('I don''t recognize the format of input_activation.  It should be a thresholded mask, clusters, or region object');
end

% process input arguments
% --------------------------------------------
doblobs = true;
dooutline = true;
doaddmontages = false;
doremove = true;
outlinecolor = [0 0 0];
splitcolor = {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]};
montagetype = 'compact';

wh = strcmp(varargin, 'noblobs');
if any(wh), doblobs = false; varargin(wh) = []; end

wh = strcmp(varargin, 'nooutline');
if any(wh), dooutline = false; varargin(wh) = []; end

wh = strcmp(varargin, 'addmontages');
if any(wh), doaddmontages = true; varargin(wh) = []; end

wh = strcmp(varargin, 'outlinecolor');
if any(wh), wh = find(wh); outlinecolor = varargin{wh(1) + 1}; end

wh = strcmp(varargin, 'splitcolor');
if any(wh), wh = find(wh); splitcolor = varargin{wh(1) + 1}; end

wh = strcmp(varargin, 'noremove');
if any(wh), doremove = false; varargin(wh) = []; end

wh = strcmp(varargin, 'full');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact2');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = false(1, length(varargin));
for i = 1:length(varargin)
    wh(i) = isa(varargin{i}, 'fmridisplay');
    if wh(i), o2 = varargin{wh}; end
end
varargin(wh) = [];

xyz = [-20 -10 -6 -2 0 2 6 10 20]';
xyz(:, 2:3) = 0;


if ~exist('o2', 'var')
    
    % set up fmridisplay
    % --------------------------------------------
    % you only need to do this once
    % then you can add montages, add and remove blobs, add and remove points (for
    % meta-analysis), etc.
    
    disp('Setting up fmridisplay objects');
    disp('This takes a lot of memory, and can hang if you have too little.');
    
    o2 = fmridisplay;
    
    % You can customize these and run them from the command line
    
    switch montagetype
        case 'full'
            o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
            o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
            
        case 'compact'
            o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
            axh = axes('Position', [0.05 0.4 .1 .5]);
            o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
        
        case 'compact2'
            o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8);
            enlarge_axes(gcf, 1.2 / 1.4);            
            axh = axes('Position', [-0.03 0.03 .2 1]);
            o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);

        otherwise error('illegal montage type. choose full or compact.')
    end
    
    wh_montages = [1 2];
    
    
else
    disp('Using existing fmridisplay object');
    
    % Other inputs will be passed into addblobs
    existingmons = length(o2.montage);
    
    if doaddmontages
        % use same o2, but add montages
        switch montagetype
            case 'full'
                o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
                o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
                
            case 'compact'
                o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6);
                axh = axes('Position', [0.05 0.4 .1 .5]);
                o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
                
            case 'compact2'
                o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8);
                enlarge_axes(gcf, 1.2 / 1.4);            
                axh = axes('Position', [-0.03 0.03 .2 1]);
                o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);

            otherwise error('illegal montage type. choose full or compact.')
        end
        
        wh_montages = existingmons + [1 2];
        
    else
        if doremove
            o2 = removeblobs(o2);
        end
        wh_montages = 1:existingmons;
  
    end
    
end

% Now we can add blobs
% --------------------------------------------

% they are added to all montages by default, but you can specify selected
% montages if you want to as well.

% it's easy to remove them as well:
% o2 = removeblobs(o2);

if doblobs
    o2 = addblobs(o2, cl, 'splitcolor', splitcolor, 'wh_montages', wh_montages, varargin{:});
end

if dooutline
    o2 = addblobs(o2, cl, 'color', outlinecolor, 'outline', 'wh_montages', wh_montages);
end


end  % function
