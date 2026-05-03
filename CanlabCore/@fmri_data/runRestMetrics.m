function results = runRestMetrics(dataObj, headerInfo, TR, varargin)
% runRestMetrics Compute resting-state metrics (ALFF, fALFF, ReHo) via the DPABI toolbox.
%
% :Usage:
% ::
%
%     results = runRestMetrics(dataObj, headerInfo, TR, ...)
%
% Computes (optionally) ALFF / fALFF and / or ReHo on a preprocessed
% fmri_data object using the DPABI toolbox
% (https://rfmri.org/DPABI). Writes the resulting maps to disk as
% NIfTI files and returns a results struct.
%
% :Inputs:
%
%   **dataObj:**
%        fmri_data object that has already been nuisance-regressed and
%        detrended.
%
%   **headerInfo:**
%        niftiinfo struct describing the target space.
%
%   **TR:**
%        Repetition time in seconds.
%
% :Optional Inputs:
%
%   **'MaskFile':**
%        Path to a mask NIfTI; the toolbox computes values only within
%        the mask. If left at default, a mask is created from the first
%        image in the input dataset. Default: 'Resampled_CanlabMask.nii'.
%
%   **'HighCutoff_ALFF':**
%        High cutoff for ALFF (Hz). Default: 0.1.
%
%   **'LowCutoff_ALFF':**
%        Low cutoff for ALFF (Hz). Default: 0.01.
%
%   **'TemporalMask':**
%        Scrubbing mask ('' = none). Default: ''.
%
%   **'ScrubbingMethod':**
%        Passed to y_alff_falff / y_reho. Default: 1.
%
%   **'NVoxel':**
%        Number of voxels in the ReHo neighbourhood. Default: 27.
%
%   **'Band_ReHo':**
%        ReHo bandpass [low high] in Hz. Default: [0.01 0.08].
%
%   **'IsNeedDetrend':**
%        ReHo detrend flag. Default: true.
%
%   **'ScrubbingTiming':**
%        Passed to y_reho. Default: [].
%
%   **'ComputeALFF':**
%        Logical, whether to compute ALFF/fALFF. Default: true.
%
%   **'ComputeReHo':**
%        Logical, whether to compute ReHo. Default: true.
%
%   **'OutDir':**
%        Folder to write .nii files into. Default: pwd.
%
%   **'ALFFFile':**
%        Filename for the ALFF output. Default: 'ALFF.nii'.
%
%   **'ReHoFile':**
%        Filename for the ReHo output. Default: 'ReHo.nii'.
%
% :Outputs:
%
%   **results:**
%        Struct with fields for each computed map and the output
%        filenames:
%
%        - .ALFFBrain / .fALFFBrain
%        - .ALFFFile
%        - .ReHoBrain
%        - .ReHoFile
%        - .Header
%
% :Examples:
% ::
%
%     % Quick outlier-detection sketch on a multi-study dataset
%     obj  = load_image_set('kragel18_alldata');
%     obj2 = rescale(obj, 'l2norm_images');
%     [est_outliers_uncorr, est_outliers_corr, outlier_tables] = ...
%         outliers(obj2, 'notimeseries');
%
%     % Example 1: Load example images, create a mask, and run rest metrics
%     fname = which('swrsub-sid001567_task-pinel_acq-s1p2_run-03_bold.nii.gz');
%     obj   = fmri_data(fname);
%     mask  = get_wh_image(obj, 1);
%
%     fname = 'mask_Example.nii';   % output file name; extension picks format
%     write(mask, 'fname', fname, 'overwrite');
%
%     % Get the TR from a sidecar JSON
%     json_struct = jsondecode(fileread(...
%         which('sub-sid001567_task-pinel_acq-s1p2_run-03_bold.json')));
%     TR = json_struct.RepetitionTime;
%
%     hdr = obj.volInfo;
%
%     results = runRestMetrics( ...
%       obj, hdr, TR, ...
%       'MaskFile'        , fname, ...
%       'NVoxel'          , 27, ...
%       'HighCutoff_ALFF' , 0.1, ...
%       'LowCutoff_ALFF'  , 0.01, ...
%       'Band_ReHo'       , [0.01 0.08], ...
%       'OutDir'          , 'C:\KeBo_Work\Restingstate_Analysistool\Result', ...
%       'ALFFFile'        , 'Alff001');
%
% :References:
%   Yan, C. G., Wang, X. D., Zuo, X. N., & Zang, Y. F. (2016). DPABI:
%   data processing & analysis for (resting-state) brain imaging.
%   Neuroinformatics, 14, 339-351.
%
% :See also:
%   - canlab_connectivity_preproc (preprocessing pipeline upstream of this)
%
% ..
%    Programmers' notes:
%    Ke Bo, 5/28/25: Added example code for running this function.
% ..

% parse inputs
p = inputParser;
p.addRequired ('dataObj'       , @(x) isa(x,'fmri_data'));
p.addRequired ('headerInfo'    , @isstruct);
p.addRequired ('TR'            , @(x) isnumeric(x)&&isscalar(x));

p.addParameter('MaskFile'       ,'Resampled_CanlabMask.nii',@ischar);
p.addParameter('HighCutoff'     ,0.1,@(x) isnumeric(x)&&isscalar(x));
p.addParameter('LowCutoff'      ,0.01,@(x) isnumeric(x)&&isscalar(x));
p.addParameter('TemporalMask'   ,'',@ischar);
p.addParameter('ScrubbingMethod',1,@(x) isnumeric(x)&&isscalar(x));

p.addParameter('NVoxel'         ,27,@(x) isnumeric(x)&&isscalar(x));
p.addParameter('Band'           ,[0.01 0.08],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('IsNeedDetrend'  ,true,@(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('ScrubbingTiming',[], @(x) isnumeric(x)||isempty(x));

p.addParameter('ComputeALFF'    ,true,@(x) islogical(x)||ismember(x,[0,1]));
p.addParameter('ComputeReHo'    ,true,@(x) islogical(x)||ismember(x,[0,1]));

p.addParameter('OutDir'         ,pwd,@ischar);
p.addParameter('ALFFFile'       ,'ALFF.nii',@ischar);
p.addParameter('ReHoFile'       ,'ReHo.nii',@ischar);

p.parse(dataObj, headerInfo, TR, varargin{:});
args = p.Results;

% ensure output dir exists
if ~exist(args.OutDir,'dir'), mkdir(args.OutDir); end

% initialize results
results = struct('Header', args.headerInfo);

% reconstruct 4-D data once
vol4D = reconstruct_image(args.dataObj);

% — ALFF/fALFF —
if args.ComputeALFF
  alffPath = fullfile(args.OutDir, args.ALFFFile);
  [ALFFBrain, fALFFBrain, ~] = y_alff_falff( ...
      vol4D, args.TR, args.HighCutoff, args.LowCutoff, ...
      args.MaskFile, alffPath, args.TemporalMask, ...
      args.ScrubbingMethod, args.headerInfo);
  results.ALFFBrain  = ALFFBrain;
  results.fALFFBrain = fALFFBrain;
  results.ALFFFile   = alffPath;
end

% — ReHo —
if args.ComputeReHo
  rehoPath = fullfile(args.OutDir, args.ReHoFile);
  [ReHoBrain, ~] = y_reho( ...
      vol4D, args.NVoxel, args.MaskFile, ...
      rehoPath, args.IsNeedDetrend, args.Band, ...
      args.TR, args.TemporalMask, args.ScrubbingMethod, ...
      args.ScrubbingTiming, args.headerInfo);
  results.ReHoBrain = ReHoBrain;
  results.ReHoFile  = rehoPath;
end
end
