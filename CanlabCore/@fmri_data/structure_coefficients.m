function sc_img_obj = structure_coefficients(weight_data, src_data_matrix)
    % structure_coefficients Compute structure coefficients (Haufe et al., 2014) for an fmri_data weight map.
    %
    % :Usage:
    % ::
    %
    %     sc_img_obj = structure_coefficients(weight_data, src_data_matrix)
    %
    % Computes structure coefficients for an fmri_data object containing
    % classifier / regression weights, given a data matrix used to fit
    % that model, following Haufe et al. (2014, NeuroImage). The
    % structure coefficient at each voxel is the correlation between the
    % voxel's data and the model's weighted predicted output, and is
    % typically more interpretable than the raw weight.
    %
    % :Inputs:
    %
    %   **weight_data:**
    %        fmri_data object whose .dat contains weight maps from a
    %        classifier or regression model (e.g., SVM, PCR).
    %
    %   **src_data_matrix:**
    %        Numeric data matrix used to fit the model. Typically raw or
    %        preprocessed voxel data, not the fmri_data object itself.
    %
    % :Outputs:
    %
    %   **sc_img_obj:**
    %        fmri_data object cloned from weight_data, with .dat
    %        replaced by the corresponding structure coefficients.
    %
    % :References:
    %   Haufe, S., Meinecke, F., Goergen, K., Daehne, S., Haynes, J.-D.,
    %   Blankertz, B., & Biessmann, F. (2014). On the interpretation of
    %   weight vectors of linear models in multivariate neuroimaging.
    %   NeuroImage, 87, 96-110.
    %
    % :See also:
    %   - structure_coefficient_map (bootstrapped voxelwise version)
    %   - get_model_encoding_map
    %
    % ..
    %    Author: Michael Sun, Ph.D.
    % ..
    
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