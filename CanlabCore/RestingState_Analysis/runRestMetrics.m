function results = runRestMetrics(dataObj, headerInfo, TR, varargin)
% runRestMetrics  Compute (optionally) ALFF/fALFF and/or ReHo on
% preprocessed data object from CanlabTool using DPABI toolbox (https://rfmri.org/DPABI)
% citation: Yan, C. G., Wang, X. D., Zuo, X. N., & Zang, Y. F. (2016). DPABI: data processing & analysis for (resting-state) brain imaging. Neuroinformatics, 14, 339-351.
%
%   results = runRestMetrics(dataObj, headerInfo, TR, ...)
%
% Inputs:
%   dataObj    – an fmri_data object (already nuisance‑regressed & detrended)
%   headerInfo – niftiinfo struct for your target space
%   TR         – repetition time (s)
%
% Name-Value Pairs:
%   'MaskFile'        – mask.nii The toolbox only compute the values within the Mask. We'll create a mask from the frist image from your image dataset (default: 'Resampled_CanlabMask.nii')
%   'HighCutoff'      – high cutoff for ALFF (Hz)   (default: 0.1)
%   'LowCutoff'       – low cutoff for ALFF (Hz)    (default: 0.01)
%   'TemporalMask'    – scrubbing mask ('' = none)  (default: '')
%   'ScrubbingMethod' – passed to y_alff_falff/y_reho (default: 1)
%   'NVoxel'          – # of voxels for ReHo        (default: 27)
%   'Band'            – ReHo bandpass [low high] Hz (default: [0.01 0.08])
%   'IsNeedDetrend'   – ReHo detrend flag           (default: true)
%   'ScrubbingTiming' – passed to y_reho (default: [])
%   'ComputeALFF'     – true/false                  (default: true)
%   'ComputeReHo'     – true/false                  (default: true)
%   'OutDir'          – folder to write .nii files  (default: pwd)
%   'ALFFFile'        – filename for ALFF output    (default: 'ALFF.nii')
%   'ReHoFile'        – filename for ReHo output    (default: 'ReHo.nii')
% :Examples:
% ::
% % -------------------------------------------------------------------------
% % Load a multi-study dataset, rescale it, and identify/plot outliers
% % Use 'notimeseries' option because this is not a time series dataset
% 
% obj = load_image_set('kragel18_alldata');
% obj2 = rescale(obj, 'l2norm_images');     % normalize heterogeneous datasets
% [est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(obj2, 'notimeseries');
%
% % -------------------------------------------------------------------------
%

%
% Output:
%   results – struct with fields for each computed map & filenames:
%     .ALFFBrain, .fALFFBrain, .ALFFFile, .ReHoBrain, .ReHoFile, .Header

% Examples and help:
% -------------------------------------------------------------------------
%
%
% Example 1: Load example images from CanlabCore and create image mask
% based on this image set
% 
% fname = which('swrsub-sid001567_task-pinel_acq-s1p2_run-03_bold.nii.gz');
% obj = fmri_data(fname);
% mask=get_wh_image(obj,1);
% 
% fname = 'mask_Example.nii';  % outname file name.  the extension (.nii, .img) determines the format.
% write(mask, 'fname', fname,'overwrite');
% 
% % Get the TR
% 
% json_struct = jsondecode(fileread(which('sub-sid001567_task-pinel_acq-s1p2_run-03_bold.json')));
% TR = json_struct.RepetitionTime;
% 
% hdr = obj.volInfo;    
% from your fmri_data
% 
% results = runRestMetrics( ...
%   obj, hdr, TR, ...
%   'MaskFile'       ,fname, ...
%   'NVoxel'         ,27, ...
%   'HighCutoff'     ,0.1, ...
%   'LowCutoff'      ,0.01, ...
%   'Band'           ,[0.01 0.08], ...
%   'OutDir'         ,'C:\KeBo_Work\Restingstate_Analysistool\Result', ...
%   'ALFFFile'       ,'Alff001');

% Programmers' notes:
% Ke Bo, 5/28/25 : Add an example code for running this function

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
