function fmri_concatenate = spm_cfg_fmri_concatenate
% SPM Configuration file for adjusting SPM.mat files after scan concatenation
%__________________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_fmri_concatenate.m 1 2022-6-3 19:26:13Z petre $


%==========================================================================
% spmmat Select SPM.mat
%==========================================================================
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {
    'Select the SPM.mat file that contains the design specification.'
    'The directory containing this file is known as the input directory.'
    }';
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

%==========================================================================
% specifications of volumes per scan
%==========================================================================
Scans          = cfg_entry;
Scans.tag      = 'scan_volumes';
Scans.name     = 'Scan Volumes';
Scans.help     = {'A 1 x n vector specifying the number of volumes in each of n scans.'
		  'Used to subdivide covariance and filter matrices into sets of matrices that are constrained to act only on certain subsets of volumes.'}';
Scans.strtype  = 'n';
Scans.num      = [1 Inf];

%--------------------------------------------------------------------------
% fmri_concatenate Model estimation
%--------------------------------------------------------------------------
fmri_concatenate          = cfg_exbranch;
fmri_concatenate.tag      = 'fmri_concatenate';
fmri_concatenate.name     = 'Concatenation Correction';
fmri_concatenate.val      = {spmmat Scans};
fmri_concatenate.help     = {
    'matlabbatch wrapper for spm_fmri_concatenate. Performs corrections to SPM.mat in the event you are using concatenated scans.'};
fmri_concatenate.prog     = @spm_run_fmri_concatenate;
fmri_concatenate.vout     = @vout_stats;
fmri_concatenate.modality = {'FMRI'};


%==========================================================================
% function dep = vout_stats(job)
%==========================================================================
function dep = vout_stats(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'SPM.mat File';
dep(1).src_output = substruct('.','spmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
