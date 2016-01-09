function new_files = remove_disdaq_vols(img_files, num_vols_per_run, num_disdaq_vols, varargin)
% :Usage:
% ::
%
%     new_files = remove_disdaq_vols(img_files, num_vols_per_run, num_disdaq_vols, ['overwrite', 0|1], ['FSLOUTPUTTYPE', fsl_output_type], ['strict', 0|1])
%
% :Inputs:
%
%   **img_files:**
%        cellstr of files to remove disdaqs from - each cell represents a run
%
%   **num_vols_per_run:**
%        vector of volume counts per run, *not* including disdaq vols
%
%   **num_disdaq_vols:**
%        constant describing how many data points to remove from the beginning of each run
%
%   **'overwrite':**
%        if set, will overwriting images in place - defaults to 0
%
%   **'FSLOUTPUTTYPE':**
%        outputfile type - defaults to 'NIFTI'
%
%   **'strict':**
%        if set, will error out unless data given to it is
%        exactly proper length - defaults to 1 - turn off ONLY with good reason
%
% :Output:
%
%   **new_files:**
%        input files, but with a 'd' prepended, unless 'overwriting' was specified
%
% NB: Not set up for 3d files yet!!!
%
% :Examples:
% ::
%
%    % for an experiment
%    num_vols_per_run = [124 140 109];  % NOT including disdaqs
%    num_disdaq_vols = 4;

    global FSLDIR;
    scn_setup();
    
    % Flags
    strict = 1;
    overwriting = 0;
    return_char = 0;
    fsl_output_type = 'NIFTI';
    
    % Check inputs
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch(varargin{i})
                    case {'overwrite' 'overwriting'}
                        overwriting = varargin{i+1};
                    case 'strict'
                        strict = varargin{i+1};
                    case {'fsl_output_type', 'FSLOUTPUTTYPE'}
                        fsl_output_type = varargin{i+1};
                end
            end
        end
    end
    
    if(ischar(img_files))
        disp('Cellstr not passed in, changing...')
        img_files = cellstr(img_files);
        return_char = 1;
    end

    % Use FSL to remove disdaqs
    new_files = cell(size(img_files));
    for i=1:length(img_files)
        current_file = img_files{i};
        num_frames = scn_num_volumes(current_file);

        if(strict && (num_frames ~= num_vols_per_run(i) + num_disdaq_vols(i)))
            error('The number of vols (%d) + disdaqs (%d) don''t add up to the length (%d) of the file (''%s'') passed in.', num_vols_per_run, num_disdaq_vols, num_frames, current_file);
        else
            if(overwriting)
                output_file = current_file;
            else
                [d f e] = fileparts(current_file);
                output_file = fullfile(d, ['d' f]);
            end

            fslroi_command = sprintf('export FSLDIR=%s && . %s/etc/fslconf/fsl.sh && FSLOUTPUTTYPE=%s && %s/bin/fslroi %s %s %d %d', ...
                FSLDIR, FSLDIR, fsl_output_type, FSLDIR, current_file, output_file, num_disdaq_vols, num_vols_per_run);
            disp(fslroi_command);
            system(fslroi_command);
        end

        new_files{i} = [output_file e];
    end

    % Convert output back to char is only a single file was passed in
    if(return_char)
        new_files = char(new_files);
    end
end
