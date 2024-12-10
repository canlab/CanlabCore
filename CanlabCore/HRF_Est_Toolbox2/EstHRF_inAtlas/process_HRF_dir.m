function [HRF_data_lvl2, HRF_data, tc_data]=process_HRF_dir(basedir, at)
    % Helper script to process an estHRF directory for second-level
    % analysis
    % Michael Sun, Ph.D.
    % - Takes a BIDS-formatted estHRF directory
    % - basedir: charstr fullpath of the directory: e.g., '//dartfs-hpc/rc/lab/C/CANlab/labdata/data/WASABI/derivatives/estHRF_NPSpos'
    % - at: atlas object
    % *Usage:
    % ::
    %    [HRF_data_lvl2, HRF_data, tc_data] = process_HRF_dir(basedir, atlas_obj})
    %    for s = 1:numel(subjects)
    %       disp(subjects{s});
    %       describeHRF(HRF_data{s});
    %       T{s}=generateHRFTable(HRF_data{s});
    %       plotHRF(HRF_data_lvl2{s}, 'FIR', my_atlas);
    %    end
    %

    % Validate input arguments
    if nargin < 2
        error('You must provide a base directory and atlas object.');
    end
    
    % Load up all the .mat files you generated with estHRF_inAtlas for each subject
    subjects = canlab_list_subjects(basedir, 'sub-*');
    HRF_data = cell(1, length(subjects));
    HRF_data_lvl2 = cell(1, length(subjects));
    tc_data = cell(1, length(subjects));
    HRF_files = cell(1, length(subjects));

    for s = 1:numel(subjects)
        HRF_files{s} = dir(fullfile(basedir, subjects{s}, '**', '*.mat'));

        % Initialize cell arrays to store the data
        HRF_data{s} = cell(1, length(HRF_files{s}));
        
        tc_data{s} = cell(1, length(HRF_files{s}));

        % Loop through each file and load the data
        for k = 1:length(HRF_files{s})
            fullpath = fullfile(HRF_files{s}(k).folder, HRF_files{s}(k).name);
            data = load(fullpath);  % Load the .mat file

            % Assuming each .mat file contains variables named 'HRF' and 'tc'
            HRF_data{s}{k} = data.HRF;
            tc_data{s}{k} = data.tc;
        end

        HRF_data_lvl2{s}=HRF_avg(HRF_data{s}, at);
    end
end