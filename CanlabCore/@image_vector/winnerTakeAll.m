function result = winnerTakeAll(input_data)
% winnerTakeAll Identify the column index of the maximum value in each voxel of an image_vector object.
%
% For each voxel (row of input_data.dat), returns the column index
% (image) holding the maximum value across images. Useful for selecting
% the most significant (largest) feature across multiple conditions or
% variables for each observation. Indices for voxels whose maximum is
% non-positive are set to 0.
%
% :Usage:
% ::
%
%    result = winnerTakeAll(input_data)
%
% :Inputs:
%
%   **input_data:**
%        An image_vector object containing the data from which the max
%        indices are to be found. If the input is a 'statistic_image',
%        it will be thresholded at p < 0.1 (uncorrected) before finding
%        the max indices.
%
% :Outputs:
%
%   **result:**
%        A structure identical to the input 'input_data', but with its
%        'dat' field replaced by a column vector of indices corresponding
%        to the maximum values in each row of the original data matrix.
%        Indices corresponding to rows where the maximum value is less
%        than or equal to 0 are set to 0.
%
% :Examples:
% ::
%
%    dat = fmri_data('results.nii'); % A 5x10 matrix of random data
%    result = winnerTakeAll(dat);
%    % result.dat will contain the indices of the max values in each row
%
% :Notes:
%
%   - If the input is a 'statistic_image', its 'dat' field is multiplied
%     by the 'sig' field after thresholding.
%   - Rows with all non-positive values in 'dat' will have a max index
%     of 0 in 'result.dat'.
%   - This function is useful for selecting the most significant
%     (largest) feature across multiple conditions or variables for each
%     observation.
%
% :See also:
%   - max
%   - threshold
%
% ..
%    Created by: Michael Sun, PhD
%    Date: 08/30/2024
%    Version: 1.0
% ..

    if isa(input_data, 'statistic_image')
        input_data = threshold(input_data, .1, 'unc');
        dat = input_data.dat .* input_data.sig;
    else
        dat = input_data.dat;
    end

    % Step 1: Initialize the result matrix with zeros
    result = input_data;

    % Step 2: Find the indices of the max values in each row
    [max_vals, max_indices] = max(dat, [], 2);

    % Step 3: Set indices where max value is less than or equal to 0 to 0
    max_indices(max_vals <= 0) = 0;

    % Step 4: Update the input_data with the indices of the max values
    result.dat = max_indices;
end
