function [dat, volInfo, cl] = iimg_threshold(image_names, varargin)
% Thresholds images to be at least thresh(1) and at most thresh(2)
%
% :Usage:
% ::
%
%     function [dat, volInfo, cl] = iimg_threshold(image_names, varargin)
%
% :Input types for image names:
%   1) String matrix of filenames
%   2) 4-D array of 3-D data volumes
%   3) voxels x images 2-D index array
%   4) image_vector object
%
% :Outputs:
%
%   **imdat:**
%        index vector of thresholded data
%
%   **volInfo:**
%        structure of info about volume
%
%   **cl:**
%        optional (slower) clusters structure from dat
%
% :Command strings:
%
%   **'imgtype':**
%        followed by 't', 'p', or 'data' (default)
%       specify type of values in image
%
%   **'threshtype':**
%        followed by 't' or 'p'
%
%   **'df':**
%        followed by degrees of freedom, for p- and t-image threshold calc
%
%   **'k':**
%        followed by extent threshold in voxels (slower)
%
%   **'abs':**
%        absolute value must be > threshold.  Use to get both pos. and neg.
%        results in two-tailed test.
%
%   **'intersect:**
%        Create intersection of images; uses abs, so + and - values
%        count as yesses for intersection
%
%   **'contrast':**
%        followed by contrasts across images
%
%   **'volInfo:**
%        followed by volInfo structure.  Necessary for extended output
%
% :Examples:
% ::
%
%    % Threshold an image set (p) to be at least zero
%    image_names =
%    /Users/tor/Documents/Tor_Documents/CurrentExperiments/Lab/Pre-appraisal/Results/rscalped_avg152T1_graymatter_smoothed.img
%    /Users/tor/Documents/Tor_Documents/CurrentExperiments/Lab/Pre-appraisal/Results/t_intercept.img
%
%    [dat, volInfo] = iimg_threshold(image_names, 'thr', 0);
%
%    % Do the same, but take the intersection and write an output image
%    [dat, volInfo] = iimg_threshold(image_names, 'thr', 3, 'outnames', 'intersct', 'masked_t.img');
%
%
%    % Threshold a t-image based on a p-value and a df
%    image_names = '/Users/tor/Documents/Tor_Documents/CurrentExperiments/Lab/Pre-appraisal/Results/t_intercept.img'
%    [dat, volInfo] = iimg_threshold(image_names, 'thr', .005, 'imgtype', 't', 'threshtype', 'p', 'df', 28, 'outnames', 'thresh_t.img');
%
%    cl = mask2clusters('masked_t.img');
%    cluster_orthviews(cl);
%
%    % The same, but threshold based on absolute value (+ and - values)
%    [dat, volInfo] = iimg_threshold(image_names, 'thr', .005, 'abs', 'imgtype', 't', 'threshtype', 'p', 'df', 28, 'outnames', 'thresh_t.img');
%
%    % Threshold a p-value image directly
%    [dat, volInfo] = iimg_threshold('X-M-Y_pvals.img', 'thr', [0 .05], 'outnames', 'X-M-Y_pvals_thresh.img');
%    [dat, volInfo, cl] = iimg_threshold(inname, 'thr', [0 .05], 'outnames', outname);
%
%    % Threshold a p-value image, with cluster sizes
%    [dat, volInfo, cl] = iimg_threshold(inname, 'thr', [0 .05], 'k', 10);
%
%    % Threshold and display a t-image using FDR, getting both positive and negative results
%    [dat, volInfo, cl] = iimg_threshold('contrast_t.img', 'imgtype', 't', 'df', 37, 'thr', .2, 'threshtype', 'fdr', 'k', 3, 'abs');
%    cluster_orthviews(cl);
%    spm_orthviews_hotcool_colormap(cat(2,cl.Z), 1.52);
%
% ..
%    Authorship and Updates:
%    Created by Tor Wager, edited by Matthew Davidson
%    Update April 2007 by TW : correct FDR-thresholding bug
% ..

    % ..
    %    Set up arguments
    % ..

    df = []; % defaults
    % maskname = deblank(image_names(1, :));
    thr = [0 Inf];                  % must be > 0, < Inf
    imgtype = 'data';               % could be data, t, p, or F
    threshtype = 'none';
    mask = [];
    k = 1;
    dointersect = 0;
    con = [];
    outnames = [];
    doabs = 0;                      % absolute value thresholding
    extended_output_flag = 2;    %*(nargout > 2); % yes if we've requested clusters to be created needed for extent also

    % inputs
    for i = 1:length(varargin)
        arg = varargin{i};
        if ischar(arg)
            switch lower(arg)
                case 'imgtype', imgtype = varargin{i+1};
                case 'threshtype', threshtype = varargin{i+1};
                case 'df', df = varargin{i+1};
                case {'thr', 'threshold'}, thr = varargin{i+1};
                case 'k', k = varargin{i+1};
                case 'intersect', dointersect = 1;
                case 'mask', mask = varargin{i+1};     
                case 'contrast', con = varargin{i+1};
                case 'outnames', outnames = varargin{i+1};
                case 'abs', doabs = 1;
                case 'volinfo', volInfo = varargin{i+1};
            end
        end
    end

    % Make upper threshold infinity if not entered
    if(length(thr) < 2)
        switch threshtype
            case 'p'
                thr = [0 thr];
            case 'fdr'
                % threshold should be a single p-value, which is converted
                % to [0 fdr_thresh] by convert_threshold, below
                % check and make sure p-value
                if thr < eps || thr > 1
                    error('Invalid threshold (thr) for FDR: must be p-value threshold.');
                end
                
            case 't'
                thr = [thr Inf];
                
            case 'none'
                % we should have a data image
                if ~strcmp(imgtype, 'data'), error('You must use imgtype ''data'' (default) if you do not enter a threshold type (''t'' or ''p'')'); end
                    
                thr = [thr Inf];
                
            otherwise
                error('Unknown threshold type (''threshtype'').');
        end
    end

    %% --------------------------------------
    % * Read image data
    % --------------------------------------
    if isa(image_names, 'statistic_image')
        tmpVolInfo = image_names.volInfo;
        image_names = replace_empty(image_names);
        dat = image_names.dat(:, 1);
    else
        [tmpVolInfo, dat] = iimg_read_img(image_names, extended_output_flag);
    end
    
    if issparse(dat), dat = full(dat); end
    clear image_names

    % use input volInfo if entered; this is so you can input an indexed image
    % instead of a filename and get extended output (cluster sizes)
    if ~exist('volInfo', 'var'), volInfo = tmpVolInfo; end
    clear tmpVolInfo

    if ~isempty(mask)
        dat = iimg_mask(mask,dat,volInfo);
    end

    if dointersect
        % If any of the input images has a zero (the designated "no value") for a voxel, then clear it out for all images
        wh = any(dat==0,2);
        dat(wh,:) = 0;
    end
    
    % --------------------------------------
    % * Convert and apply threshold
    % --------------------------------------
    thr = convert_threshold(imgtype, threshtype, thr, df, dat);

    if doabs
        whzero = abs(dat) < thr(1) | abs(dat) > thr(2);  % remove these values
    else
        whzero = dat < thr(1) | dat > thr(2);
    end
    dat(whzero) = 0;

    %% --------------------------------------
    % * Apply size threshold
    % --------------------------------------
    [dat, nvox] = iimg_cluster_extent(dat, volInfo, k);


    %% --------------------------------------
    % special operations
    % --------------------------------------
    if ~isempty(con)
        % Apply contrast(s) across images
        dat = apply_contrast(dat, con);
    end

    %% --------------------------------------
    % * Write out image names, if asked
    % --------------------------------------
    if ~isempty(outnames)
        iimg_write_images(dat, volInfo, outnames);
    end

    %% --------------------------------------
    % * clusters output
    % --------------------------------------
    if extended_output_flag

        % could enter u and k here, but don't need to because it's already
        % thresholded; also, neg. elements will be removed in some cases if
        % threshold is re-applied here.
        
        %cl = iimg_indx2clusters(dat, volInfo, thr, k);
        cl = iimg_indx2clusters(dat, volInfo);
    end
end




% Sub-functions

function thr = convert_threshold(imgtype, threshtype, thr, df, dat)
    % handle FDR
    switch imgtype
        case 'data'
            if(strcmp(threshtype, 'fdr'))
                error('FDR thresholding has no meaning for data images. Try using a t- or a p-image');
            end
        case 't'
            switch threshtype
                case 't'
                case 'p'
                    % convert p-value threshold to t-value threshold
                    if ~exist('df', 'var') || isempty(df)
                        error('You must enter df to use p-value based thresholding with a t img.');
                    end
                    thr(~isinf(thr)) = tinv(1-thr(~isinf(thr)), df);
                case 'fdr'
                    % convert t-values in image to p-values for FDR
                    % thresholding
                    if ~exist('df', 'var') || isempty(df)
                        error('You must enter df to use p-value based thresholding with a t img.');
                    end
                    pdat = 1 - tcdf(dat(dat ~= 0 & ~isnan(dat)), df); % zero is a special invalid value

                    fdr_threshold = FDR(pdat, thr); % get p-value threshold
                    if isempty(fdr_threshold), fdr_threshold = 0; end
                    
                    fprintf('Computed FDR threshold. For p < %.3f corrected, p-threshold is %.7f, ', thr, fdr_threshold);
                    thr = [tinv(1-fdr_threshold, df) Inf];
                    
                    fprintf('t-threshold is %.2f\n', thr(1));
                    
                otherwise
                    error('Threshold type must be t or p values for t-images')
            end

        case 'p'
            switch threshtype
                case 't'
                    if ~exist('df', 'var') || isempty(df)
                        error('You must enter df to use t-value based thresholding with a p img.');
                    end
                    thr = 1 - tcdf(thr, df);
                case 'p' % do nothing
                case 'fdr'
                    
                    dat = dat(dat ~= 0 & ~isnan(dat)); % zero is a special invalid value
                    
                    fdr_threshold = FDR(dat, thr);   
                    if isempty(fdr_threshold), fdr_threshold = 0; end
                    
                    %thr(isinf(thr)) = fdr_threshold;
                    fprintf('Computed FDR threshold. For p < %.3f corrected, p-threshold is %.7f, ', thr, fdr_threshold);
                    
                    thr = [0 fdr_threshold];  % from zero to new thresh
                otherwise
                    error('Threshold type must be t or p values for p-images')
            end
    end
end


function condat = apply_contrast(dat, con)
    if size(con, 2) ~= size(dat, 2)
        error('Contrast matrix must have as many columns as images.');
    end

    condat(:, i) = dat * con';
end



