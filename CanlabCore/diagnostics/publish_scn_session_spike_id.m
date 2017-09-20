function publish_scn_session_spike_id(inputimgs, SUBJDATA,outformat)
% This function is a wrapper function to call scn_session_spike_id in
% 'multi-session' mode, using input data across the runs for a single
% subject.  It runs the program, and generates both a yaml-format text file
% for uploading into the CANlab database, and an html file with all the
% results and images for that subject embedded.
%
% :Inputs:
%
%   **inputimgs:**
%        is a cell array of images (4-D) for each run in a separate cell.
%        the cell array outside this function has to be named 'inputimgs'
%
%   **SUBJDATA:**
%        Input fields of SUBJDATA define the experiment name, subject name,
%        and directories for saving both QC images + yaml and HTML
%        - study. a string 
%        - subject. another string with subject ID
%        - html_save_dir can be one directory for all subjects if format 'pdf'
%        - subject_dir. a new subdirectory qc_images will be created within
%        subject_dir. images and .yaml files will be written there.
%
%
%   **ouformat**
%       file format for the report as available in 'publish'. can be
%       'html' [default],'pdf', 'latex', etc.
%
% :Examples:
% ::
%
%    SUBJDATA.study = 'NSF';
%    SUBJDATA.subject = subjects{i};
%    SUBJDATA.html_save_dir = fullfile(output_basedir, 'html_output');
%    SUBJDATA.subject_dir = fullfile(output_basedir, 'SubjectData', 'denoised_canlab', SUBJDATA.subject);
%
% ..
%    This has been tested with SPM8 and Matlab2010a only.
%
%    Tor Wager, Aug 2010
% ..

% Programmer's Notes
%
% 8/25/2017 - Stephan Geuter
% added option to use other fileformats like PDF
%

%% check input
% default report format is html. choose pdf or other options from publish.m
if nargin<3, outformat = 'html'; end


%% Initialize yaml file for database integration
if ~exist(SUBJDATA.subject_dir, 'dir'), mkdir(SUBJDATA.subject_dir); end
cd(SUBJDATA.subject_dir)

yamlfilename = fullfile(SUBJDATA.subject_dir, 'qc_results.yaml');

SUBJDATA.unique_id = [SUBJDATA.study '_' SUBJDATA.subject];

struct2yaml(yamlfilename, SUBJDATA, 'new', 'replace');

%% housekeeping: get rid of any old .png qc images
% these will mess it up, as it checks for existing files in order to
% require minimal input to scn_session_spike_id

qcdir = fullfile(SUBJDATA.subject_dir, 'qc_images');
if exist(qcdir, 'dir')
    disp('Removing old qc png images');
    
    disp(['!rm ' qcdir filesep '*png'])
    eval(['!rm ' qcdir filesep '*png'])
end

%% Create .m file script to publish that will create HTML file

htmlname = ['qc_scn_session_spike_id_' SUBJDATA.subject];

if ~exist(SUBJDATA.html_save_dir, 'dir'), mkdir(SUBJDATA.html_save_dir); end

diaryfilename = fullfile(SUBJDATA.html_save_dir, [htmlname '.m']);
if exist(diaryfilename, 'file'), eval(['!rm ' diaryfilename]); end

% These commands, written in the file called diaryfilename,
% will be evaluated in the m-file, and image and text
% results published in the HTML file called htmlname.
% ------------------------------------------------------------------

diary(diaryfilename);

disp(['%% SCN session spike ID, subject ' SUBJDATA.subject]);
fprintf('%% *Experiment:* %s\n%% \n', SUBJDATA.study); % start new formatted block if two spaces
fprintf('%% *QC results stored in:* %s\n%% \n', SUBJDATA.subject_dir);
fprintf('%% spike ID, quality control step 1.  writes to ''qc_results.yaml''\n');
%%% NEED COMMAND THAT LOADS INPUTIMGS FROM DISK -- OR PASSES IN AS INPUT?
disp('scn_session_spike_id(inputimgs);')
disp('snapnow');

diary off

%% 
addpath(SUBJDATA.html_save_dir)
p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 500, ...
    'outputDir', SUBJDATA.html_save_dir, 'showCode', false,'format',outformat);
publish(diaryfilename, p);

fprintf('\npublished to\n%s\n\n',diaryfilename);

if strcmp(outformat,'html')
    htmlfilename = fullfile(SUBJDATA.html_save_dir, [htmlname '.html']);
    web(htmlfilename);
end
    

% clean up the original diary; no longer needed
eval(['!rm ' diaryfilename]);



    
    

    
