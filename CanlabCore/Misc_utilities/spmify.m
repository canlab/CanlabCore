function gKWYs=spmify(obj, SPM)
   % spmify: Transform input data using SPM (Statistical Parametric Mapping) processing
    %
    % This function takes input data (in various forms) and an SPM structure,
    % and outputs SPM gKWY time series data as a double matrix. It conforms the
    % data to the expected number of runs based on the SPM model and applies
    % necessary SPM filters and transformations.
    %
    % Inputs:
    %   obj - Input data. Can be a file path (or list of file paths), an 'fmri_data' object,
    %         or a data matrix (single or double precision).
    %   SPM - A structure containing SPM design matrix and related data.
    %
    % Outputs:
    %   gKWYs - Processed time series data in the form of a double matrix, adjusted
    %           according to the SPM structure provided.
    %
    % Notes:
    %   - The function checks for the data type of 'obj' and processes it accordingly.
    %   - If more than one run is expected based on the SPM design matrix, the function
    %     will error out if only one object is passed unless run information can be
    %     reliably extracted from the input data.
    %   - This function assumes that run numbers in BIDS-formatted filenames (if used)
    %     conform to the concatenation order of the SPM design matrix.
    %
    % Example:
    %   SPM = load('SPM.mat');
    %   data = 'path/to/fmri_data.nii';
    %   processedData = spmify(data, SPM);
    %
    % See also: fmri_data, spm_filter
    % Michael Sun, Ph.D. 02/09/2024


    K=SPM.xX.K;
    runs_expected=numel(K);

    % Check to see if obj conforms to runs_expected, otherwise wisely
    % select what to extract from SPM from SPM metadata

    if strcmpi(class(obj), 'fmri_data') || ischar(obj) || isstring(obj)
        if runs_expected > 1
            % Extract run number and use that to match that with the
            % appropriate g, K, and W parameters in SPM
            % CAUTION: This may not work if run numbers assigned do not
            % conform to the concatenation number of the SPM design matrix.
            % bidsInfo=extractBIDSinfo(obj.image_names);
            % run=bidsInfo('run');
            error('More than 1 run expected, but only 1 object passed in.');
        else
            switch class(obj)
                case 'fmri_data'
                    y{1}= obj.dat;
                case 'char'
                    obj=fmri_data(obj);
                    y{1}= obj.dat;
                case 'string'
                    obj=fmri_data(char(obj));
                    y{1}= obj.dat;
            end
        end
    elseif iscell(obj)
        for i=1:numel(obj)
            switch class(obj{i})
                case 'fmri_data'
                    y{i}= obj{i}.dat;
                case 'char'
                    obj{i}=fmri_data(obj{i});
                    y{i}= obj{i}.dat;
                case 'string'
                    obj{i}=fmri_data(char(obj{i}));
                    y{i}= obj{i}.dat;
                case 'single'
                    y{i}=obj{i};
                case 'double'
                    y{i}=obj{i};
            end
        end
    end

    if strcmpi(class(obj), 'single') || strcmpi(class(obj), 'double') 
        if size(obj, 3)~=runs_expected
            error('Dimensions of matrix data does not match number of runs expected.')
        end
        for i=1:numel(size(obj, 3))
            y{i}=obj(:,:,i);
        end
    end

    % Sort your data in a cellarray
    for i=1:runs_expected
        W{i}=SPM.xX.W(SPM.xX.K(i).row, SPM.xX.K(i).row);
        KWY{i}=spm_filter(K(i),W{i}*double(y{i}'));
    end

    gKWYs = cell(1, 4);
    
    offset = 0;
    for j = 1:length(KWY)
        for i = 1:size(KWY{j}, 1)
            gKWYs{j}(:, i) = KWY{j}(i, :) * SPM.xGX.gSF(i + offset);
        end
        offset = offset + size(KWY{j}, 1);
    end
    

end