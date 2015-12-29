function MASKSTATS = canlab_maskstats(msks,imgs,varargin)
% Produces comparison of pattern mask and images
% e.g., look at NPS pattern expression in set of beta images
%
% :Usage:
% ::
%
%    MASKSTATS = canlab_maskstats(maskfiles,imgfiles,[options])

% :Inputs:
%
%   **maskfiles:**
%      string or cellstring of mask filenames or fmri_data object

%   **imgfiles:**
%      string or cellstring of image filenames or fmri_data object
%
% :Optional Inputs:
%
%   **ts:**
%        timeseries treatment: each string in imgfiles is assumed to be
%        a 4D file. Data will be returned with one column per time series and
%        one volume per row. If not all timeseries are same length, all will
%        be NaN-padded to the length of the longest timeseries.
%
%        Note: does not work with imgfiles input as fmri_data object
%
%   **dir2cell:**
%        will sort stats into cells based on directory containing the imgfile
%        they belong to such that each cell contains one directory's worth of
%        stats which is a vector with a value for each imgfile.
%        :Examples:
%        ::
%
%             # Input: a list of single trial betas for a set of subjects
%             b = filenames('sub*/*heat_trials*.img');
%             ms = canlab_maskstats('nps',b,'dot_product','dir2cell');
%
%        Output: includes the set of cells that are input to a mediation
%        analysis.
%
%   **keepzeros:**
%        don't remove zeros from imgfiles before taking measurements
%
%   **keepzerosmask:**
%        don't remove zeros from maskfiles before taking measurements
%
%   **single:**
%        leave data as single (DEFAULT: convert to double)
%
%   **trinarize:**
%        trinarize maskfile (set values larger than 0 to 1 and values less than
%        zero to -1)
%
%   **noreshape:**
%        don't attempt to reshape results according to imgfiles array
%
%   **nobin:**
%        don't binarize mask before extracting mean or std
%
% :Note: ts, dir2cell, and noreshape are mutually exclusive options
%
% **Built-In Masks**
% The following strings can be given as the maskfile argument to
% call up built-in mask files:
%
%   **nps:**
%        weights_NSF_grouppred_cvpcr.img
%
%   **nps_thresh:**
%        weights_NSF_grouppred_cvpcr_FDR05.img
%
%   **nps_thresh_smooth:**
%        weights_NSF_grouppred_cvpcr_FDR05_smoothed_fwhm05.img
%
% :Measure Options:
%
%   **all:**
%        add: mean, dot_product, centered_dot_product,
%             cosine_similarity, and correlation
%
%   **mean:** (DEFAULT)
%      apply binarized mask to images and return means
%      mean(img .* abs(bin(mask)))
%
%   **std:**
%        apply binarized mask to images and return standard deviations
%        std(img .* abs(bin(mask)))
%
%   **dot_product:**
%        dot(mask, img)
%
%   **cosine_similarity:**
%        dot(mask, img) / (norm(mask) * norm(img))
%
%   **correlation:**
%        corr(mask, img)
%
%   **centered_dot_product:**
%        dot(mask-mean(mask), img-mean(img))
%
%
% :Details:
%   - imgfiles are spatially resampled to maskfiles
%   - voxels with zeros in maskfile are removed
%   - in-mask voxels with zeros in imgfiles will generate warnings
%
% ..
%    AUTHOR: Luka Ruzic 2013
% ..


OP = {}; % struct to contain desired "operations" (e.g., dot product, mean, etc)
%% parse arguments
ALLOPS = false; % adds all implemented operations
DO_TRINARIZE = false; % signed binarize option
DO_ZERO2NAN_DATA = true; % treat zeros as non-data in data input
DO_ZERO2NAN_MASK = true; % treat zeros as non-data in mask input
DO_RESHAPE = true; % attempt to reshape output data according to shape of input data
DO_BIN = true; % binarize mask
DO_TS = false; % timeseries mode
DO_CELLS = false; % cell array output mode

i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}   
            case {'binmean' 'mean' 'norm' ...
                  'dot_product' 'centered_dot_product' 'cosine_similarity' 'correlation'}
                OP{end+1} = varargin{i}; %#ok    
            case 'all'
                ALLOPS = true;
            case 'keepzeros'
                DO_ZERO2NAN_DATA = false;
            case 'keepzerosmask'
                DO_ZERO2NAN_MASK = false;
            case 'trinarize'
                DO_TRINARIZE = true;
            case 'nobin'
                DO_BIN = false;
            case 'noreshape'
                DO_RESHAPE = false;
                DO_CELLS = false;
                DO_TS = false;
            case 'ts'                
                DO_RESHAPE = true;
                DO_TS = true;
                DO_CELLS = false;
            case 'dir2cell'
                DO_RESHAPE = true;
                DO_CELLS = true;                
                DO_TS = false;
            otherwise
                error('Unrecognized argument %s',varargin{i})
        end
    elseif iscellstr(varargin{i})
        for j=1:numel(varargin{i})
            switch varargin{i}{j}
                case {'binmean' 'mean' 'std' 'nonbinmean' 'nonbinstd' 'norm' ...
                        'dot_product' 'centered_dot_product' 'cosine_similarity' 'correlation'}
                    OP{end+1} = varargin{i}{j}; %#ok
                case 'all'
                    ALLOPS = true;
                otherwise
                    error('Unrecognized argument %s',varargin{i}{j})
            end
        end
    else
        disp(varargin{i})
        error('Above argument unrecognized')
    end
    i=i+1;
end

if ALLOPS, OP = [OP 'mean' 'std' 'dot_product' 'centered_dot_product' 'cosine_similarity' 'correlation']; end
    
if isempty(OP), OP = {'mean'}; end


%% error checking
if isempty(msks), error('Must provide one or more maskfiles'); end
if ischar(msks), msks = cellstr(msks); end
if ~iscellstr(msks) && ~isa(msks,'fmri_data'), error('maskfiles must be a string, cell array of strings, or fmri_data object');end

if isempty(imgs), error('Must provide one or more imgfiles'); end
if ischar(imgs), imgs = cellstr(imgs); end
if ~iscellstr(imgs) && ~isa(imgs,'fmri_data'), error('imgfiles must be string, cell array of strings, or fmri_data object'); end

MASKSTATS(numel(msks)) = struct; % initializing


%% prepare mask files
if ~isa(msks,'fmri_data')
    for m = 1:numel(msks)
        switch msks{m}
            case 'nps'
                msks{m} = which('weights_NSF_grouppred_cvpcr.img');
            case 'nps_thresh'
                msks{m} = which('weights_NSF_grouppred_cvpcr_FDR05.img');
            case 'nps_thresh_smooth'
                msks{m} = which('weights_NSF_grouppred_cvpcr_FDR05_smoothed_fwhm05.img');
            otherwise
                if ~exist(msks{m},'file')
                    error('No such file: %s',msks{m})
                end
        end
    end
end


%% prepare data
% fprintf('LOADING IMAGE FILES\n')
if isa(imgs,'fmri_data')
    imgdata = imgs;
    imgdata = imgdata.replace_empty;
    if DO_TS
        warning('''ts'' mode is incompatible with fmri_data input: turning off'); %#ok
        DO_TS = false;
    end    
    imgfiles = cellstr(imgdata.fullpath);
else
    imgfiles = imgs;
    evalc('imgdata = fmri_data(imgfiles);');       
end

imgdata.dat = double(imgdata.dat); % convert singles to doubles
if DO_ZERO2NAN_DATA, imgdata.dat(imgdata.dat==0) = NaN; end  % remove zeros if desired
    
if DO_TS
    for i = 1:numel(imgfiles)
        ni = nifti(imgfiles{i});
        tsdim(i) = ni.dat.dim(4); %#ok
    end    
end

if DO_CELLS    
%     [ignore1 ignore2 x] = unique(regexprep(imgfiles,'/[^/]*$','')); %#ok
    [ignore1 ignore2 x] = unique(regexprep(cellstr(imgdata.fullpath),'/[^/]*$','')); %#ok
    reshaperows = histc(x,1:max(x));    
end


%% produce comparison
if isa(msks,'fmri_data')
    allmaskdata = msks;
    allmaskdata = allmaskdata.replace_empty;
    nmasks = size(allmaskdata.dat,2);
else
    nmasks = numel(msks);
end
for m = 1:nmasks
    clear imgdat maskdat imgdatcent maskdatcent imgdatnorm maskdatnorm
    
    % save filenames
    % imgs (modify cell array according to mode)
    if DO_CELLS
        for i = 1:numel(reshaperows)
            st = sum(reshaperows(1:i-1))+1;
            MASKSTATS(m).imgfiles{i} = imgfiles(st:st+reshaperows(i)-1);
        end
    elseif DO_TS        
        MASKSTATS(m).imgfiles = imgfiles';
    else
        MASKSTATS(m).imgfiles = imgfiles;
    end    
    
    % load and prep mask data

    if isa(msks,'fmri_data')
        maskdata = allmaskdata;        
        maskdata.dat(:,[1:size(allmaskdata.dat,2)]~=m) = [];
        MASKSTATS(m).maskfile = deblank(allmaskdata.fullpath(m,:));
    else
        evalc('maskdata = fmri_data(msks{m});');
        MASKSTATS(m).maskfile = msks{m};
    end
    fprintf('GETTING DATA FOR MASK: %s\n',MASKSTATS(m).maskfile);
    maskdata.dat = double(maskdata.dat); % convert singles to doubles
    %if DO_ZERO2NAN_MASK, maskdata.dat(maskdata.dat == 0) = NaN; end % SHOULD THIS HAPPEN HERE?   remove zeros
    maskdata = maskdata.resample_space(imgdata); % resample
    if DO_ZERO2NAN_MASK, maskdata.dat(maskdata.dat == 0) = NaN; end % remove zeros
    
    % select in-mask voxels
    wh_inmask = ~isnan(maskdata.dat);
    maskdat = maskdata.dat(wh_inmask);
    imgdat = imgdata.dat(wh_inmask,:);    
    
    % trinarize
    if DO_TRINARIZE
        maskdat(maskdat<0) = -1;
        maskdat(maskdat>0) = 1;
    end
    
    % warn about in-mask zeros
    if DO_ZERO2NAN_DATA
        inmaskzeros = isnan(imgdat);
        if any(inmaskzeros)
            fprintf('WARNING! There are zeros in your data within this mask!\n')
            %fprintf('Consider carefully whether this impairs your ability to compare data points.\n')
            %fprintf('  These are being ignored as non-data but your results will be impacted!\n')
            for i = find(sum(inmaskzeros))
                fprintf(' %7d in-mask zeros in %s\n',sum(inmaskzeros(:,i)),deblank(imgdata.fullpath(i,:)));
            end
        end
    end
    
    % perform operation(s)
    for i = 1:numel(OP)
        switch OP{i}
            case 'mean'
                if DO_BIN
                    MASKSTATS(m).stats.mean = nanmean(imgdat)';
                else
                    MASKSTATS(m).stats.nonbinmean = nanmean(bsxfun(@times,imgdat,maskdat))';
                end
                
            case 'std'
                if DO_BIN
                    MASKSTATS(m).stats.std = nanstd(imgdat)';
                else                    
                    MASKSTATS(m).stats.nonbinstd = nanstd(bsxfun(@times,imgdat,maskdat))';
                end
                
            case 'norm'
                MASKSTATS(m).stats.norm = (nansum(imgdat .^ 2) .^ .5)';
                
            case 'dot_product'                
                MASKSTATS(m).stats.dot_product = nansum(bsxfun(@times,imgdat,maskdat))';
                
            case 'centered_dot_product'                
                imgdatcent = bsxfun(@minus,imgdat,nanmean(imgdat));
                maskdatcent = bsxfun(@minus,maskdat,nanmean(maskdat));
                MASKSTATS(m).stats.centered_dot_product = nansum(bsxfun(@times,imgdatcent,maskdatcent))';
                
            case 'cosine_similarity'
                imgdatnorm = nansum(imgdat .^ 2) .^ .5;
                maskdatnorm = nansum(maskdat .^ 2) .^ .5;
                MASKSTATS(m).stats.cosine_similarity = (nansum(bsxfun(@times,imgdat,maskdat)) ./ (imgdatnorm .* maskdatnorm))';                              
                
            case 'correlation'
                MASKSTATS(m).stats.correlation = corr(imgdat,maskdat,'rows','pairwise');
                
            otherwise
                error('Comparison type (%s) unrecognized',OP{i})
        end
        
        % reshape according to mode
        if DO_RESHAPE
            if DO_TS
                if numel(unique(tsdim)) == 1
                    tslen = tsdim(1);
                else
                    % not all timeseries are same length
                    % assume all start at same time and pad with NaNs
                    tslen = max(tsdim);
                    for j = 1:numel(tsdim)
                        cutind = ((j-1)*tslen) + tsdim(j);
                        MASKSTATS(m).stats.(OP{i}) = [... 
                            MASKSTATS(m).stats.(OP{i})(1:cutind); ...
                            nan(tslen-tsdim(j),1); ...
                            MASKSTATS(m).stats.(OP{i})(cutind+1:end)];
                    end 
                    % alternative coding option: vector -> cells -> NaN-pad the cells -> cell2mat %
                end
                MASKSTATS(m).stats.(OP{i}) = reshape(MASKSTATS(m).stats.(OP{i}),[tslen numel(imgfiles)]);
            elseif DO_CELLS                
                MASKSTATS(m).stats.(OP{i}) = mat2cell(MASKSTATS(m).stats.(OP{i}),reshaperows,1)';
            else                
                try MASKSTATS(m).stats.(OP{i}) = reshape(MASKSTATS(m).stats.(OP{i}),size(imgfiles)); end %#ok
            end
        end
    end    
end

end
