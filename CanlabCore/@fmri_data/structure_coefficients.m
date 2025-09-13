function sc_img_obj = structure_coefficients(weight_data, src_data_matrix)
    % COMPUTE_STRUCTURE_COEFFICIENTS computes structure coefficients for
    % an fMRI data object using a given data matrix following Haufe et al., 2014.
    %
    % Inputs:
    %   weight_data - fMRI data object containing weight data the .dat from
    %   classification e.g., SVM or PCR
    %   src_data_matrix - data matrix to correlate, usually raw or
    %   preprocessed data. Not the fmri_data object itself
    %
    % Output:
    %   sc_img_obj - fMRI data object with structure coefficients in sc_img_obj.dat
    %
    %   Michael Sun, Ph.D.
    
    % Calculate the correlation matrix
    corr_matrix = corr(src_data_matrix');
    
    % Initialize the structure coefficients matrix with the same size as img_obj.dat
    structure_coefficients = zeros(size(weight_data.dat));
    
    % Compute structure coefficients by multiplying the correlation matrix with
    % each column of weight_data.dat
    for i = 1:size(weight_data.dat, 2)
        structure_coefficients(:, i) = corr_matrix * weight_data.dat(:, i);
    end
    
    % Create the output fMRI data object with structure coefficients
    sc_img_obj = weight_data;  % Copy the original object structure
    sc_img_obj.dat = structure_coefficients;  % Replace the data with structure coefficients
end