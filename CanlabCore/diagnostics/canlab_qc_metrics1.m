function [QC, dat] = canlab_qc_metrics1(epi_names, mask_name, varargin)
% Calculate quality control metrics on a 4-D EPI Analyze or Nifti file
%
% Standard CANlab usage is to enter a single 4-D ravol* for one run, and
% the brain mask implicit_mask.img created in canlab_preproc.m
%
% :Inputs:
%
%   **epi_names:**
%        Names of (usually one) 4-D EPI file, in cell or string, full path
%
%   **mask_name:**
%        Name of brain mask image, string, full path
%              IF EMPTY: Uses implict masking (better) and calculates
%              ghost/signal outside mask
%
% :Optional Inputs:
%
%   **noplot:**
%        skip plots
%
%   **noverbose:**
%        skip output printout to screen
%
%   **printfile:**
%        followed by name of file to print output to, full path
%
%   **noheader:**
%        suppress printing of header (var names) in output
%
% :Missing values and basic info:
%
%   **num_images:**
%        number of images
%
%   **missing_vox:**
%        Voxels in mask with NaN values or zero values at every time point
%
%   **missing_images:**
%        Images with NaN values or zero values at every voxel
%
%
%   **missing_values:**
%        NaN or zero values in valid images / voxels. Zeros could
%        be interpreted as values of zero in analysis, causing artifacts in results
%        if these are actually invalid values.
%
% Missing voxels will often be ignored in analyses in most software, but
% Missing images/values could cause problems
%
% :Basic signal to noise:
%
%
%   **perc-mean-ghost:**
%        mean signal outside the mask / mean total signal
%
%
%   **mean_snr:**
%        mean Cohen's d (signal/noise, SNR) across time (temporal SNR) within the mask.
%    Mean divided by standard deviation across time at each voxel, averaged.  Higher is better.
%
%
%   **snr_inhomogeneity:**
%        standard deviation of SNR within the mask. Lower is better.
%
%
%   **snr_inhomogeneity95:**
%        95% confidence range for SNR within the mask. Lower is better.
%
%
%   **rms_successive_diffs:**
%        Essentially a high-pass filtered version of SNR,
%        expressed as a fraction of the overall mean and averaged across voxels. Lower is better.
%
%
%   **rms_successive_diffs_inhomogeneity:**
%        standard deviation of the above across voxels. Lower is better.
%
% :Left-right asymmetry:
%
%
%   **signal_rms_asymmetry:**
%        Voxel-wise left/right root mean square asymmetry in mean signal across time, expressed as
%        a fraction of the mean SNR.  Reflects both gross inhomogeneity and noise.  Lower is better.
%
%
%   **signal_hemispheric_asymmetry:**
%        Root mean square difference between left and
%        right hemispheres, expressed as a fraction of the grand mean signal across time.  Reflects
%        gross inhomogeneity.  Lower is better.
%
%
%   **snr_rms_asymmetry:**
%        Voxel-wise left/right root mean square asymmetry in SNR, expressed as
%        a fraction of the mean SNR.  Reflects both gross inhomogeneity and noise.  Lower is better.
%
%
%   **snr_hemispheric_asymmetry:**
%        Root mean square difference between left and
%        right hemispheres, expressed as a fraction of the mean SNR.  Reflects
%        gross inhomogeneity.  Lower is better.
%
% :Examples:
% ::
%
%    %SETUP
%    mydir{1} = '/Users/tor/Documents/Tor_Documents/Coursework_and_Teaching/psyc7215/Data/UM_Face_House/060518mw/Functional/Preprocessed/run_01';
%    wcard = 'rarun*img';
%    epi_names = filenames(fullfile(mydir{1}, wcard), 'absolute');
%
%    maskdir = '/Users/tor/Documents/Tor_Documents/Coursework_and_Teaching/psyc7215/Data/UM_Face_House/060518mw';
%    mask = 'implicit_mask.img';
%    mask_name = fullfile(maskdir, maskname);
%
%    %RUN
%    QC = canlab_qc_metrics1(epi_names, mask_name);
%
%    QC = canlab_qc_metrics1(epi_names, mask_name, 'noplot', 'printfile', 'test_qc.txt');
%    QC = canlab_qc_metrics1(epi_names, mask_name, 'noplot', 'printfile', 'test_qc.txt', 'noheader');
%

doplot = 1;
verbose = 1;
QC = [];
printfile = [];
noheader = 0;

for i = 1:length(varargin)
    
    if ischar(varargin{i})
        
        switch varargin{i}
            % functional commands
            case 'noplot', doplot = 0;
            case 'noverbose', doverbose = 0;
            case 'noheader', noheader = 1;
            case 'printfile', printfile = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


if iscell(epi_names), epi_names = strvcat(epi_names{:}); end

% idstr = epi_names';
% idstr = idstr(:);
% if length(idstr) > 100, idstr = idstr(1:100); end
idstr = epi_names;
if size(idstr,2) > 50, idstr = spm_file(idstr(1,:),'short25'); end



QC = struct('filenames', epi_names(:), 'idstr', idstr);
%all_outputs = {'filenames'};
all_outputs = {};

% Get implicit mask, if no mask is entered
if nargin < 3 || isempty(mask_name)
 
    mask_name = fmri_mask_image(epi_names); % full mask with all valid voxels
    
%     maskobj = fmri_mask_image(epi_names, 'implicit');
%     indx = maskobj.volInfo.image_indx;
%     indx(maskobj.volInfo.wh_inmask(maskobj.removed_voxels)) = false;

end

% Load data
% ------------------------------------------
dat = fmri_data(epi_names, mask_name);

% Implicit mask
% ------------------------------------------
if isa(mask_name, 'image_vector')
    
    fprintf('Calculating implicit mask: ');
    m = mean(dat);
    [dummy, dummy, nvox, is_inmask] = fmri_mask_thresh_canlab(m); %to use dip-based mask: fmri_mask_thresh_canlab(m, [],'dip');
    fprintf('%3.0f voxels in, %3.0f voxels out of mask.\n', nvox, sum(~is_inmask));
    
    out_of_mask_dat = dat.dat(~is_inmask, :);
    
    dat = remove_empty(dat, ~is_inmask, []);
    
    perc_mean_ghost = mean(out_of_mask_dat(:)) ./ mean(dat.dat(:));
    
end


% do this here; we will use in plot
snr = mean(dat.dat') ./ std(dat.dat');

% Plot data
% ------------------------------------------
if doplot
    plot(dat); 

    h = findobj('Tag', 'fmri data matrix');
    figure(h);
    subplot(2, 3, 6);
    [h, x] = hist(snr, 100);
    h = h ./ sum(h);
    plot(x, h, 'k', 'LineWidth', 3);
    
    title('SNR distribution (PDF)');
end


% Missing values
% ------------------------------------------

outnames = {'num_images', 'vox_vol', 'missing_vox', 'missing_images', 'missing_values'};
if exist('perc_mean_ghost', 'var'), outnames = [outnames 'perc_mean_ghost']; end
    
all_outputs = [all_outputs outnames];

dat = remove_empty(dat);

vox_dims = diag(dat.volInfo.mat(1:3, 1:3));
vox_vol = prod(abs(vox_dims));
num_images = size(dat.removed_images, 1);
missing_vox = sum(dat.removed_voxels);
missing_images = sum(dat.removed_images);

missing_values = sum(all(isnan(dat.dat) | dat.dat == 0, 2));

for i = 1:length(outnames)
    
    eval(['myvar = ' outnames{i} ';']);
    QC.(outnames{i}) = myvar;
    
end

% SNR
% ------------------------------------------
outnames = {'mean_snr', 'mean_snr_per_mm', 'snr_inhomogeneity', 'snr_inhomogeneity95' 'rms_successive_diffs' 'rms_successive_diffs_inhomogeneity'};
all_outputs = [all_outputs outnames];

%snr = mean(dat.dat') ./ std(dat.dat');
mean_snr = mean(snr);
mean_snr_per_mm = mean_snr ./ vox_vol;

snr_inhomogeneity = std(snr);
snr_inhomogeneity95 = prctile(snr, 95) - prctile(snr, 5);

rmssd = ( diff(dat.dat') .^ 2 ./ size(dat.dat, 2) ) .^ .5; % rms at each voxel
rms_successive_diffs = mean(rmssd);

rms_successive_diffs_inhomogeneity = std(rmssd);

for i = 1:length(outnames)
    
    eval(['myvar = ' outnames{i} ';']);
    QC.(outnames{i}) = myvar;
    
end

% Asymmetry calculation
% ------------------------------------------
outnames = {'signal_rms_asymmetry', 'signal_hemispheric_asymmetry', 'snr_rms_asymmetry', 'snr_hemispheric_asymmetry'};
all_outputs = [all_outputs outnames];

dat = replace_empty(dat);
m = mean(dat.dat');
snr = m ./ std(dat.dat');

[signal_rms_asymmetry, signal_hemispheric_asymmetry] = calculate_asymmetry(dat, m);

[snr_rms_asymmetry, snr_hemispheric_asymmetry] = calculate_asymmetry(dat, snr);

for i = 1:length(outnames)
    
    eval(['myvar = ' outnames{i} ';']);
    QC.(outnames{i}) = myvar;
    
end


% Print
% ------------------------------------------
for i = 1:length(all_outputs)
    out(1, i) = mean(QC.(all_outputs{i}));
end

if ~isempty(printfile)
    diary(printfile);
end

if verbose || printfile
    if noheader
        % no header row
        print_matrix(out, [], {QC.idstr});
    else
        print_matrix(out, all_outputs, {QC.idstr})
    end
    
end

if ~isempty(printfile)
    diary off
end

end % main function




function [rms_asymmetry, hemispheric_asymmetry] = calculate_asymmetry(dat, snr)

xyz = round(voxel2mm(dat.volInfo.xyzlist', dat.volInfo.mat)');
isleft = xyz(:, 1) < 0;
isright = xyz(:, 1) > 0;
xyz(:, 1) = abs(xyz(:, 1));
[u1, indx1] = unique(xyz, 'first', 'rows');
[u2, indx2] = unique(xyz, 'last', 'rows');

if any(any(u1 - u2)), disp('Bug in asymmetry algorithm'); end
wh_unique = indx1 ~= indx2; % these have two distinct voxels with the same mm coordinates after rounding and abs x

lrdat = [snr(indx1(wh_unique)); snr(indx2(wh_unique))];

rms_asymmetry = sqrt( nansum(diff(lrdat) .^ 2) ./ size(lrdat, 2) ) ./ nanmean(lrdat(:));

tmpsnr = snr;
tmpsnr(~isleft) = NaN;
leftsnr = tmpsnr([indx1(wh_unique); indx2(wh_unique)]);
leftsnr = nanmean(leftsnr);

tmpsnr = snr;
tmpsnr(~isright) = NaN;
rightsnr = tmpsnr([indx1(wh_unique); indx2(wh_unique)]);
rightsnr = nanmean(rightsnr);

hemispheric_asymmetry = sqrt((leftsnr - rightsnr)^2 ./ 2) ./ sum([leftsnr rightsnr]);

end


