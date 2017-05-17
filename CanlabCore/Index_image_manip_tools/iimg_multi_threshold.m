function [cl, dat, cl_extent, dat_extent] = iimg_multi_threshold(inname, varargin)
% :Usage:
% ::
%
%     [cl, dat, cl_extent, dat_extent] = iimg_multi_threshold(inname, varargin)
%
% Multi-threshold application and orthview display of blobs
%
% :Inputs:
%
%   **inname:**
%        Either the name of an image file or a data vector
%
%        if data, enter volInfo structure separately as below.
%
% :Command strings:
%
%   **'prune':**
%        consider only contiguous blobs for which at least 1 voxel meets
%        the most stringent threshold
%
%   **'pruneseed':**
%        followed by a vectorized thresholded activation image
%
%        Only those blobs overlapping with at least one 'seed' data voxel
%        are saved
%
%        Also: the first color (highest threshold) in the output images is assigned to the seed
%
%   **'add':**
%        add new blobs to existing orthviews
%
%   **'p':**
%        if image is a p-value image (and thresholds are p-values)
%
%
%   **'thresh':**
%        followed by a vector of n thresholds, most to least stringent
%
%   **'size':**
%        followed by a vector of size thresholds
%
%   **'volInfo': followed by volInfo struct; needed if inname is data rather
% than filenames
%
%   **'overlay':**
%        followed by overlay image (anatomical) filename
%
%   **'transseed':**
%        transparent seed
%
%   **'hideseed':**
%        do not show seed regions on plot
%
%   **'pos':**
%        positive results only
%
%   **'neg':**
%        negative results only
%
%   **'add':**
%        add to existing multi-threshold plot
%
% Needs (development): Legend, nice handling of colors, input
% colors, color maps specified with command words (like 'red')
%
% :Examples:
% ::
%
%    inname = 'Activation_proportion.img';
%    [cl, dat] = iimg_multi_threshold(inname, 'prune', 'thresh', [.1 .5 .3], 'size', [1 5 10]);
%
%    cl2 = iimg_multi_threshold(Pimg(1, :), 'thresh', [.001 .01 .05], 'size', [3 5 10], 'p');
%    cl2 = iimg_multi_threshold(Pimg(1, :), 'thresh', [.001 .01 .05], 'size', [3 3 3], 'p', 'prune');
%
% from act + corr results (see robust_results_act_plus_corr)
% First prepare 'seed' regions that overlap with correlated regions, then
% use multi_threshold to see full extent of clusters
% ::
%
%    [dattype, seeddat] = iimg_check_indx(res.act.overlapdat, res.volInfo, 'full');
%    [cl, dat] = iimg_multi_threshold('rob_p_0001.img', 'p', 'thresh', [.001 .01 .05], 'size', [1 1 1], 'pruneseed', seeddat)
%
% Display an F-map from robust regression on a customized mean anatomical,
% with pruning.
% ::
%
%    cl = iimg_multi_threshold('rob_pmap_full.img', 'thresh', [.001 .01 .05], 'size', [1 1 1], 'p', 'prune', 'overlay', EXPT.overlay);
%
% Display regions at 3 thresholds with an input data 'seed' vector
% ::
%
%    [cl, dat] = iimg_multi_threshold(pvals, 'p', 'thresh', [.001 .01 .05], 'size', [1 1 1], 'pruneseed', p_sig', 'volInfo', R.volInfo);
%
% As above, but ONLY POSITIVE RESULTS
% ::
%
%    [cl, datout] = iimg_multi_threshold(pimg1, 'thresh', [.001 .01 .05], 'size', [1 1 1], 'p', 'pruneseed', dat, 'overlay', EXPT.overlay, 'colors', colors, 'transseed', 'pos', 'rob_tmap_0001.img');
%
% Complete example for showing positive and negative blobs:
% ::
%
%    [clpos] = iimg_multi_threshold('slope_p.img', 'p', 'thresh', [.005 .01 .05], 'size', [1 1 1], 'prune', 'overlay', anat, 'pos', 'slope_t.img');
%    [clneg] = iimg_multi_threshold('slope_p.img', 'p', 'thresh', [.005 .01 .05], 'size', [1 1 1], 'prune', 'overlay', anat, 'neg', 'slope_t.img', 'add', 'colors', {[0 0 1] [0 .5 1] [.4 .5 1]});
%
% ..
%    tor wager, July 1, 2006
%
%    Programmers' notes:
%    tor: changed to SPM8 overlay, April 2011
%    tor: fixed bug in table printing (cluster_extent) that was not using sign
%    info to calculate extent, returning unsigned extent maps identical across
%    positive and negative results in some cases.
% ..


% ..
%    Setup and defaults
% ..

cl = [];
cl_extent = [];
dat = [];
dat_extent = [];
signdat = [];

% spm 99 overlay = which('scalped_single_subj_T1.img');
overlay = which('SPM8_colin27T1_seg.img');  % spm8 seg cleaned up
if isempty(overlay)
    overlay = which('spm2_single_subj_T1_scalped.img');
end

mask = [];
thresh = [.10 .05 .04];
szthresh = [1 10 10];
colors = [1 1 0; 1 .5 0; 1 .3 .3; 0 1 0; 0 0 1];
dotransseed = 0;
dohideseed = 0;
add2existing = 0;
dopruneclusters = 0;
ispimg = 0;
posnegstr = 'both';
signimg = [];
wh_handle = 1;
fdr_thresh = .05;       % default FDR threshold

for i = 1:length(varargin)
    if ischar(varargin{i}) && ~isempty(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'prune', dopruneclusters = 1;
            case 'add', add2existing = 1;
            case 'p', ispimg = 1;
                
                % functional commands
            case 'pruneseed', dopruneclusters = 1; pruneseed = varargin{i+1};
            case 'noprune', dopruneclusters = 0;
                
            case 'thresh', thresh = varargin{i+1};
            case 'size', szthresh = varargin{i+1};
            case 'volInfo', volInfo = varargin{i+1};
            case 'overlay', overlay = varargin{i+1}; varargin{i+1} = [];
            case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
            case 'colors', colors = varargin{i+1};
            case 'transseed', dotransseed = 1;
            case 'hideseed', dohideseed = 1;
            case 'pos', posnegstr = 'pos'; signimg = varargin{i+1}; varargin{i+1} = []; % image containing signs of vectors
            case 'neg', posnegstr = 'neg'; signimg = varargin{i+1}; varargin{i+1} = [];
                
            case 'fdrthresh', fdr_thresh = varargin{i + 1};
                
            case 'handle', wh_handle = varargin{i+1};
                
            case 'noadd' % do nothing
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ischar(inname) && ~isempty(inname)
    % deal with 4-D files: use only first volume
    inname = check_valid_imagename(inname, 1);
    
    inname = expand_4d_filenames(inname);
    if size(inname, 1) > 1
        fprintf('Using only first volume in 4-D image file: \n%s\n', inname(1,:));
        inname = inname(1, :);
    end
    
elseif isempty(inname)
    disp('Input image name is empty!')
    return
end

if exist('signimg', 'var') && ~isempty(signimg) && ischar(signimg)
    signimg = check_valid_imagename(signimg, 1);
    
    % deal with 4-D files: use only first volume
    signimg = expand_4d_filenames(signimg);
    if size(signimg, 1) > 1
        fprintf('Using only first volume in 4-D image file: \n%s\n', signimg(1,:));
        signimg = signimg(1, :);
    end
end

num_thresholds = length(thresh);

fprintf(1, '\niimg_multi_threshold viewer\n=====================================================\n');
fprintf(1, 'Showing positive or negative: %s\n', posnegstr);
fprintf(1, 'Overlay is: %s\n', overlay);
if isempty(mask)
    fprintf(1, 'No mask specified in iimg_multi_threshold\n');
elseif ischar(mask)
    fprintf(1, 'Mask is: %s\n', mask);
else
    fprintf(1, 'Mask is entered in iimg_multi_threshold as a vector of numerical values\n');
end

ynstr = {'No' 'Yes'};
fprintf(1, 'Entered p-value image? : %s\n', ynstr{ispimg+1});
fprintf(1, 'Height thresholds:')
fprintf(1, '%3.4f ', thresh);
fprintf(1, '\n');
fprintf(1, 'Extent thresholds: ');
fprintf(1, '%3.0f ', szthresh);
fprintf(1, '\n');
ynstr = {'No' 'Yes'};
fprintf(1, 'Show only contiguous with seed regions: %s\n ', ynstr{dopruneclusters+1});
fprintf(1, '\n');



% --------------------------------------
% Load image with signed values, if we have one
% --------------------------------------
if ~isempty(signimg)
    
    if ischar(signimg)
        [dummy, signdat] = iimg_read_img(signimg);
        
    elseif isa(signimg, 'image_vector')
        signdat = signimg.dat(:, 1);
    
    elseif isa(signimg, 'double')
        signdat = signimg;    
        
    else
        error('Enter string name or image_vector object for sign image.');
    end
    
    switch posnegstr
        case 'pos'
            signdat = signdat > 0;
        case 'neg'
            signdat = signdat < 0;
        otherwise
            error('Sign string should be pos or neg.');
    end
end

% --------------------------------------
% Apply mask now, to simplify things later
% Also apply mask from sign image, if we have one
% Now mask has all relevant masking data
% --------------------------------------
% Tor edited 7/2014

if ~isempty(mask) && ischar(mask)
    mask = iimg_read_img(mask);
end

if ~isempty(signdat)
    if isempty(mask), mask = true(size(signdat)); end
    
    if any(size(signdat) ~= size(mask))
        error('Mask and sign data image do not appear to be the same dims.  Check inputs.');
    end
    mask = mask .* signdat;
end

% --------------------------------------
% Threshold images - FIRST threshold only
% --------------------------------------
if ispimg
    current_thresh = [0 thresh(1)];
    lowest_thresh = [0 thresh(end)];    % for pruning based on extent
    
    disp('Warning: p-values in cl.Z will not give valid spm_max subclusters.')
    disp('log(1/p) saved in output cl.Z field.');
    
    imgtype = 'p';
    if(thresh(1) == Inf)
        threshtype = 'fdr';
        current_thresh = fdr_thresh; % default FDR or user-specified
    else
        threshtype = 'p';
    end
else
    if exist('signimg', 'var') && ~isempty(signimg) && ischar(signimg) && strcmp(posnegstr, 'neg')
        % reverse to go from -Inf to cutoff for negative responses only
        current_thresh = [-Inf thresh(1)];
        lowest_thresh = [-Inf thresh(end)];
    else
        current_thresh = [thresh(1) Inf];
        lowest_thresh = [thresh(end) Inf];
    end
    imgtype = 'data';
    threshtype = 'none';
end

% Actually load the image and threshold
% --------------------------------------
if isa(inname, 'statistic_image')
    dat = threshold(inname, current_thresh(2), 'unc', 'k', szthresh(1));
    volInfo = dat.volInfo;
    dat = replace_empty(dat);
    dat = dat.dat(:, 1);    
       
elseif exist('volInfo', 'var')
    % use input volInfo if entered; this is so you can input an indexed image
    % instead of a filename and get extended output (cluster sizes)
    dat = iimg_threshold(inname, 'thr', current_thresh, 'k', szthresh(1), 'imgtype', imgtype, 'threshtype', threshtype, 'mask', mask, 'volInfo', volInfo);

elseif ischar(inname)
    % This option if inputs are filenames
    [dat, volInfo] = iimg_threshold(inname, 'thr', current_thresh, 'k', szthresh(1), 'imgtype', imgtype, 'threshtype', threshtype, 'mask', mask);
else
    error('Enter statistic_image object or string file name for input image.');
end


% restrict to only pos or neg values if posneg option is specified
if ~isempty(signimg)
    dat = dat .* signdat;  % all data, which will be used as pruneseed
    
    if exist('pruneseed', 'var')
        pruneseed = pruneseed .* signdat;
    end
end


if dopruneclusters
    if exist('pruneseed', 'var')
        % find vox at Lowest threshold that are contiguous with seed
        % this defines 'extent'.  Then, later, get any voxels that are
        % contiguous with the voxels in the extent region at higher thresholds.
        % this way, voxels at an intermediate thresh that are contiguous with
        % the extent region but not with the actual seed are included.
        % --------------------------------------
        % tor edit: 8/1/09.  add mask.  when masking, extent region can
        % cover masked areas, leaving 'orphan' voxels that do not obey
        % the extent threshold
        % tor edit: 7/22/14. mask should now incorporate signdat info as well.

        dat_extent = get_extent_regions(inname, pruneseed, lowest_thresh, szthresh(end), volInfo, mask);
    else
        % use highest threshold as seed
        dat_extent = get_extent_regions(inname, dat(:,1), lowest_thresh, szthresh(end), volInfo, mask);
    end
    
    cl_extent = iimg_indx2clusters(dat_extent, volInfo);
    
    if sum( abs(dat_extent) > eps*10 ) > 50000
        disp('Warning! Too many voxels at lowest threshold to cluster.  Cluster extent-based pruning will not work correctly.');
    end
    
else
    % aug 2010 fix: if only one threshold, clpos_extent == clpos, so
    % mediation tables will print correctly
    dat_extent = get_extent_regions_noprune(inname, lowest_thresh, szthresh(end), volInfo, mask);
    cl_extent = iimg_indx2clusters(dat_extent, volInfo);
end




% --------------------------------------
% Get other thresholds
% --------------------------------------

dat = [dat zeros(size(dat, 1), length(thresh)-1)];

if dopruneclusters  % exist('pruneseed', 'var') % Tor changed, 8/1/09
    % save only blobs that show some activation in extent region (lowest thresh) around seed input image
    dat(:,1) = iimg_cluster_prune(dat(:,1), dat_extent, volInfo);
end

for i = 2:num_thresholds
    
    if ispimg
        current_thresh = [0 thresh(i)];
        
        imgtype = 'p';
        if(thresh(i) == Inf)
            threshtype = 'fdr';
        else
            threshtype = 'p';
        end
    else
        current_thresh = [thresh(i) Inf];
        imgtype = 'data';
        threshtype = 'none';
    end
    
    
    if thresh(i) >= 0
        dat(:,i) = iimg_threshold(inname, 'thr', current_thresh, 'k', szthresh(i), 'imgtype', imgtype, 'threshtype', threshtype, 'mask', mask, 'volInfo', volInfo);
    else
        % negative threshold; invalid for p-maps
        dat(:,i) = iimg_threshold(inname, 'thr', [-Inf thresh(i)], 'k', szthresh(i), 'imgtype', imgtype, 'threshtype', threshtype, 'mask', mask, 'volInfo', volInfo);
    end
    
    if dopruneclusters
        % --------------------------------------
        % save only blobs that show some activation
        % WITHIN extent region ; Tor, 8/1/09
        % --------------------------------------
        dat(:, i) = dat(:, i) .* abs(dat_extent) > 10*eps;
        
        %             if exist('pruneseed', 'var')
        %                 dat(:,i) = iimg_cluster_prune(dat(:,i), dat_extent, volInfo);
        %             else
        %                 dat(:,i) = iimg_cluster_prune(dat(:,i), dat_extent, volInfo);
        %             end
    end
end


% --------------------------------------
% Load image with signed values, if we have one
% Apply sign to dat
% --------------------------------------
if ~isempty(signimg)
    % restrict to only pos or neg values if posneg option is specified
    dat = dat .* repmat(signdat, 1, num_thresholds);
end




% --------------------------------------
% Append 'seed', if entered, as "highest threshold" image vector
% So it gets shown and assigned the first color
% --------------------------------------
if exist('pruneseed', 'var')
    if dotransseed
        % don't include seed in things to plot w/solid colors
        % this also means don't remove voxels that are in seed regions from
        % lower thresholds
        clseed = iimg_indx2clusters(pruneseed, volInfo);
    else
        dat = [pruneseed dat];
        num_thresholds = num_thresholds+1;
    end
end

% --------------------------------------
% Report significant voxels
% --------------------------------------

% Get rid of voxels that appear in earlier ("higher") maps
mysum = cumsum(dat, 2);
mysum(:, end) = [];
mysum = [zeros(size(mysum, 1), 1) mysum];
mysum = mysum > 0;
dat(mysum) = 0;

% --------------------------------------
% Report significant voxels
% --------------------------------------
sigvox = sum(dat>0);

% make clusters
cl = cell(1, num_thresholds);
for i = 1:num_thresholds
    
    %voldata = iimg_reconstruct_3dvol(dat(:,i), volInfo);
    %cl{i} = mask2clusters(voldata, volInfo.mat);
    
    cl{i} = iimg_indx2clusters(dat(:,i), volInfo);
end


% --------------------------------------
% define colors
% --------------------------------------
%colors = [linspace(1, 0, 3)' linspace(1, .3, 3)' linspace(0, .3, 3)'];


% --------------------------------------
% image clusters on brain (orthviews)
% --------------------------------------

if add2existing, addstr = 'add'; else addstr = 'noadd'; end

if exist('pruneseed', 'var')
    if dotransseed && ~dohideseed
        cluster_orthviews(clseed, {colors(end, :)}, 'overlay', overlay, 'trans', addstr, 'handle', wh_handle, 'skipempty');
        addstr = 'add';
    end
end

for i = 1:num_thresholds
    if i == 1 && dohideseed && ~dotransseed
        % dat(:,1) is seed, but we have asked not to show it
        % do nothing
    else
        if iscell(colors)
            thiscolor = colors{i};
        else
            thiscolor = colors(i,:);
        end
        
        cluster_orthviews(cl{i}, {thiscolor}, 'overlay', overlay, 'solid', addstr, 'handle', wh_handle, 'skipempty');
        if ~isempty(cl{i}), addstr = 'add'; end % after first one, add additional to existing
        
        if ispimg
            % replace p-values with log(1/p) in Z field.  Higher is thus more
            % sig.  for compatibility with spm_max
            for j = 1:length(cl{i})
                cl{i}(j).Z_descrip = 'Log(1/p)';
                cl{i}(j).Z = log(1./cl{i}(j).Z);
            end
        end
        
    end
end


% Title on figure
fh = findobj('Tag', 'Graphics');
ch = get(fh, 'Children');
if isempty(ch)
    disp('Cannot find SPM graphics window. It should have a tag of ''Graphics''.');
else
    for i= 1:length(ch), mytype = get(ch(i), 'Type'); wh(i)=strcmp(mytype, 'axes'); end
    if ~exist('wh','var')
        disp('Found SPM graphics window, but cannot find any axes.');
    else
        axish = ch(wh);
        %axish = sort(axish);
        axes(axish(1));
    end
end

if exist('pruneseed', 'var') && ~dotransseed
    seedstr = 'seed ';
else
    seedstr = ' ';
end
if(any(isinf(thresh)))
    thrstr = repmat(' %3.6f', 1, length(thresh)-1);
    str = sprintf(['%s thr: %s FDR' thrstr], posnegstr, seedstr, thresh(~isinf(thresh)));
else
    thrstr = repmat(' %3.6f', 1, length(thresh));
    str = sprintf(['%s thr: %s' thrstr], posnegstr, seedstr, thresh);
end

title(str, 'FontSize', 18);
end




function dat_extent = get_extent_regions(inname, pruneseed, lowest_thresh, szthresh, volInfo, mask)
if lowest_thresh >= 0
    dat = iimg_threshold(inname, 'thr', lowest_thresh, 'k', szthresh, 'volInfo', volInfo, 'mask', mask);
else
    % negative threshold; invalid for p-maps
    dat = iimg_threshold(inname, 'thr', [-Inf lowest_thresh(1)], 'k', szthresh, 'volInfo', volInfo, 'mask', mask);    
end

dat_extent = iimg_cluster_prune(dat, pruneseed, volInfo);
end


function dat_extent = get_extent_regions_noprune(inname, lowest_thresh, szthresh, volInfo, mask)

if lowest_thresh >= 0
    dat = iimg_threshold(inname, 'thr', lowest_thresh, 'k', szthresh, 'volInfo', volInfo, 'mask', mask);
else
    % negative threshold; invalid for p-maps
    dat = iimg_threshold(inname, 'thr', [-Inf lowest_thresh(1)], 'k', szthresh, 'volInfo', volInfo, 'mask', mask);
end

dat_extent = dat;
end
