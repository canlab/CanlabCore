function stats_image_obj = threshold(stats_image_obj, input_threshold, thresh_type, varargin)
%
% stats_image_obj = threshold(stats_image_obj, pvalthreshold or other thresh, thresh_type, ['k', extent_thresh])
%
% Inputs:
%
% input_threshold: [pvalthreshold or other thresh]
% A numeric value corresponding to the threshold desired.
% Either a p-value or a range of raw values, depending on the threshold
% type.
%
% thresh_type: Threshold type
% can be one of:
% 'fdr' : FDR-correct based on p-values already stored in image .p field
% 'fwe' : FWE-correct; not implemented
% 'bfr' : Bonferroni correction (FWE).
% 'unc' : Uncorrected p-value threshold: p-value, e.g., .05 or .001
%
% 'raw-between' : threshold raw image values; save those > input_threshold(1) and < input_threshold(2)
% 'raw-outside' : threshold raw image values; save those < input_threshold(1) or > input_threshold(2)
%
% 'k', followed by cluster extent in voxels: extent-based thresholding of
% any of the above
%
% Tor Wager, Dec 2010
% 
% 7/19/2013: added 'mask' and 'bfr'(bonferroni) option by Wani Woo. 
%   With the 'mask' option, you can define a space for the multiple comparison 
%   correction. 
% e.g. dat = threshold(dat, 0.001, 'unc', 'k', 35, 'mask', which('scalped_avg152T1_graymatter_smoothed.img'));
%      dat = threshold(dat, 0.001, 'unc', 'k', 35, 'mask', maskobj);

k = 1;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'k'
                k = varargin{i + 1};
                
                if isempty(stats_image_obj.volInfo)
                    error('You must add a volInfo structure to the statistic image object to do extent-based thresholding');
                end
                
            case 'mask' % 7/19/2013: added by Wani Woo. 
                
                maskinput = varargin{i + 1}; varargin{i + 1} = 0;
                switch class(maskinput)
                    case 'char' % string file name
                        maskobj = fmri_mask_image(maskinput);
                        
                    case 'fmri_mask_image'
                        maskobj = maskinput;
                        
                    otherwise
                        error('region class constructor: unknown mask input type.')
                end
                clear maskinput
                
                if length(maskobj.dat) ~= length(stats_image_obj.dat)
                    if sum(sum(maskobj.volInfo.mat ~= stats_image_obj.volInfo.mat))
                        mask_mapdat = scn_map_image(maskobj, stats_image_obj);
                        maskdat = logical(mask_mapdat(stats_image_obj.volInfo.wh_inmask));
                    else
                        error('VolInfo of mask and data are same, but data lengths are different.');
                    end 
                else
                    if sum(sum(maskobj.volInfo.mat ~= stats_image_obj.volInfo.mat))
                        error('VolInfo of mask and data are diffferent, but data lengths are same.');
                    else
                        maskdat = logical(maskobj.dat);
                    end
                end
                
                stats_image_obj = replace_empty(stats_image_obj);
                stats_image_obj.p = stats_image_obj.p(maskdat);
                stats_image_obj.dat = stats_image_obj.dat(maskdat);
                if ~isempty(stats_image_obj.sig)
                    stats_image_obj.sig = stats_image_obj.sig(maskdat);
                end
                stats_image_obj.removed_voxels = ~maskdat;
                
            otherwise
                error('Illegal argument keyword.');
        end
    end
end


n = size(stats_image_obj.dat, 2); % images are columns here

for i = 1:n
    
    % Determine threshold
    % ---------------------------------------------------
    
    switch lower(thresh_type)
        case 'fdr'
            % fix zero P-vals, which can occur with some nonparam methods
            % or large t-values
            wh = stats_image_obj.p == 0 & abs(stats_image_obj.dat) > 1000*eps;
            stats_image_obj.p(wh) = 1000*eps;
            
            thrval = FDR(stats_image_obj.p(:, i), input_threshold);
            if isempty(thrval)
                thresh(i) = 0;
            else
                thresh(i) = thrval;
            end
            fprintf('Image %3.0f FDR q < %3.3f threshold is %3.6f\n', i, input_threshold, thresh(i));
            stats_image_obj.threshold(i) = thresh(i);
            
        case 'fwe'
            error('not implemented yet.')
        
        case {'bfr', 'bonferroni'} % 7/19/2013: added by Wani Woo. 
            thresh(i) = input_threshold/size(stats_image_obj.p,1);
            stats_image_obj.threshold(i) = thresh(i);
            
        case {'unc', 'uncorrected'}
            thresh(i) = input_threshold;
            stats_image_obj.threshold(i) = thresh(i);
            
        case {'raw-between', 'raw-outside'}
            thresh = input_threshold;
            
        otherwise
            error('Unknown threshold type. Enter fdr, fwe, or uncorrected, or raw-between or raw-outside')
    end
    
    %     % Re-set .dat depending on image type (needs to be updated for other
    %     % image types)
    %     % ---------------------------------------------------
    %     switch
    
    switch lower(thresh_type)
        
        case {'fdr', 'fwe', 'unc', 'uncorrected', 'bfr', 'bonferroni'}
            stats_image_obj.sig(:, i) = stats_image_obj.p(:, i) < thresh(i);
            
            % Apply size threshold
            % --------------------------------------
            stats_image_obj = replace_empty(stats_image_obj);
            if k > 1
                stats_image_obj.sig(:, i) = logical(iimg_cluster_extent(double(stats_image_obj.sig(:, i)), stats_image_obj.volInfo, k));
            end
            
            fprintf('\nImage %3.0f\nPositive effect: %3.0f voxels, min p-value: %3.8f\n', i, ...
                sum(stats_image_obj.dat(:, i) > 0 & stats_image_obj.sig(:, i)), ...
                min(stats_image_obj.p(stats_image_obj.dat(:, i) > 0, i)));
            
            fprintf('Negative effect: %3.0f voxels, min p-value: %3.8f\n', ...
                sum(stats_image_obj.dat(:, i) < 0 & stats_image_obj.sig(:, i)), ...
                min(stats_image_obj.p(stats_image_obj.dat(:, i) < 0, i)));
            
        case 'raw-between'
            stats_image_obj.sig(:, i) = stats_image_obj.dat(:, i) > thresh(1) & stats_image_obj.dat(:, i) < thresh(2);
            
            % Apply size threshold
            % --------------------------------------
            if k > 1
                stats_image_obj.sig(:, i) = logical(iimg_cluster_extent(double(stats_image_obj.sig(:, i)), stats_image_obj.volInfo, k));
            end
            
            fprintf('Keeping vals between %3.3f and %3.3f: %3.0f voxels in .sig\n', thresh(1), thresh(2), sum(stats_image_obj.sig(:, i)));
            
        case 'raw-outside'
            stats_image_obj.sig(:, i) = stats_image_obj.dat(:, i) < thresh(1) | stats_image_obj.dat(:, i) > thresh(2);
            
            % Apply size threshold
            % --------------------------------------
            if k > 1
                stats_image_obj.sig(:, i) = logical(iimg_cluster_extent(double(stats_image_obj.sig(:, i)), stats_image_obj.volInfo, k));
            end
            
            fprintf('Keeping vals outside of %3.3f to %3.3f: %3.0f voxels in .sig\n', thresh(1), thresh(2), sum(stats_image_obj.sig(:, i)));
            
        otherwise
            error('Unknown threshold type. Enter fdr, fwe, or uncorrected')
    end
    
    % Replace .dat field with zero for non-sig voxels, so that we can
    % create regions, etc.  .t field and .p field, etc., still have
    % original info.
    % But doing this here breaks the ability to re-threshold, so now done
    % in region.m
    %stats_image_obj.dat(~stats_image_obj.sig(:, i), i) = 0;
    
    
end  % image loop

if isempty(stats_image_obj.volInfo)
    disp('Warning! volInfo not defined for stats image object.  regions cannot be formed correctly from this image.')
else
    
    
    % Update volInfo, for region-defining purposes.
    % Note that this will only work for the FIRST image in the set
    
    stats_image_obj = reparse_contiguous(stats_image_obj);
    
    % i = 1;
%     wh = logical(stats_image_obj.sig(:, i));
%     n = sum(wh);
%     
%     if n < 50000
%         stats_image_obj.volInfo(1).cluster = zeros(size(stats_image_obj.volInfo(1).cluster));
%         stats_image_obj.volInfo(1).cluster(wh) = spm_clusters(stats_image_obj.volInfo(1).xyzlist(wh, :)')';
%     else
%         stats_image_obj.volInfo(1).cluster(wh) = ones(n, 1);
%     end
end

end % function



