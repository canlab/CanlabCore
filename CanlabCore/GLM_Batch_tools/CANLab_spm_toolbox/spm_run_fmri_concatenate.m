function out = spm_run_fmri_est(job)
% Estimate parameters of a specified model
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_fmri_est.m 7354 2018-06-22 10:44:22Z guillaume $


%-Load SPM.mat file
%--------------------------------------------------------------------------
load(job.spmmat{1},'SPM');
if ~exist('SPM','var')
    error('The MAT-file does not contain an SPM variable.');
end
out.spmmat = job.spmmat;

%-Move to the directory where the SPM.mat file is
%--------------------------------------------------------------------------
original_dir = pwd;
cd(fileparts(job.spmmat{:}));

spm_fmri_concatenate(job.spmmat{1}, job.scan_volumes)

%out.spmvar = SPM;
cd(original_dir);
fprintf('Done\n')
