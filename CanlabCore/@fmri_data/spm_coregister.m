function out_obj = spm_coregister(obj, varargin)
% Use SPM to coregister an image to a standard MNI template or another image
%
% - writes files to disk using SPM, in same directory as loaded from
% - need reference images on disk; will not work (yet) with objects only
% - also realigns if 4D
%
% out_obj = spm_coregister(obj, ['template', template_image_name])
%

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

template = '/Users/torwager/Documents/GitHub/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/icbm152_2009_symm_enhanced_for_underlay.img,1';

doplot = false;
doverbose = true;

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'template'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Make sure template is char array with filename
if iscell(template)
    template = template{1};
elseif isa(template, 'image_vector')
    template = template.fullpath;
end

job.ref = {template};

job.source = cellstr(obj.fullpath); % {'/Users/torwager/Dropbox (Dartmouth College)/COURSES/Courses_Dartmouth/2020_2_Winter_fMRI_Class/Shared_resources_for_students/PSYC60_class_MRI_datalad/1071_Psych60_Pinel_Movie/sub-sid001567/anat/sub-sid001567_acq-MPRAGE_T1w.nii,1'};
job.other = {''};

% If source (to move) is a 4-D file, coreg first vol only
job.source = cellstr(expand_4d_filenames(job.source));

if size(job.source, 1) > 1
    source_is_4d = true;
    
    % If we use only the first image, SPM will apparently move all other
    % images along with it. If we include all images as source images, it
    % will separately coreg each. Adding to "other" does not seem
    % necessary, as it will apply coreg to all frames of 4D image for each
    % other frame, potentially. 
    % OR MAYBE NOT...need all images or spaces are messed up
    %
    % for future: coreg to mean if 4d, don't reslice, then realign and
    % reslice
    
    job.other = job.source(2:end, :);
    
    job.source = job.source(1);
    
else
        source_is_4d = false;
end

% (check)
% spm_check_registration(strvcat(job.source, job.other(1:5)));

    
% Set options

job.eoptions = struct( ...
'cost_fun', 'nmi', ... % 'tol' [1×12 double],
'fwhm', [4 4]); % [7 7]

job.roptions = struct( ...
'interp', 4, ...
'wrap', [0 0 0], ...
'mask', 0, ...
'prefix', 'r');

% x  = spm_coreg(char(job.ref), char(job.source), job.eoptions);

out = spm_run_coreg(job);

if source_is_4d
    % Strip final image frame specification (,***)
    [mypath, myfile, ext] = fileparts(out.rfiles{1});
    out.rfiles = fullfile(mypath, [myfile ext(1:4)]);
    
    % realign
    Vout = spm_realign(out.rfiles);
    spm_reslice(Vout)
    
    % get realigned file
    out.rfiles = fullfile(mypath, ['r' myfile ext(1:4)]);
end
 
% Create object
out_obj = fmri_data(out.rfiles, 'noverbose');

end
