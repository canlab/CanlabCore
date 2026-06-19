function obj = check_properties(obj)
% check_properties Check and fill in empty properties for a statistic_image object.
%
% Fill in empty fields (.p, .ste, .sig) of a statistic_image object if
% needed, expanding stored data to the full in-mask voxel grid using
% volInfo.wh_inmask.
%
% ***under construction, do not use yet***
%
% :Usage:
% ::
%
%     obj = check_properties(obj)
%
% :Inputs:
%
%   **obj:**
%        A statistic_image object whose .p, .ste, and/or .sig fields may
%        be empty. The function references obj.volInfo.nvox and
%        obj.volInfo.wh_inmask to expand fields onto the full image grid.
%
% :Outputs:
%
%   **obj:**
%        The input object with empty fields populated to the full voxel
%        grid (default values: p = 1, ste = Inf, sig = 0 for voxels
%        outside the mask).
%
% :See also:
%   - statistic_image
%   - replace_empty
%   - validateattributes
%
% ..
%    2018 July : Created by Tor Wager
% ..

    obj_out = replace_empty(obj_out);
    k = size(obj_out.dat, 2);
    
    for i = 1:k
        % this may break if nvox (total in image) is different for 2
        % images...
        
        % Wani: in some cases, obj could have empty p, ste, and sig
        if ~isempty(obj.p)
            p = ones(obj.volInfo.nvox, k);
            p(obj.volInfo.wh_inmask, i) = obj.p(:, i);
        end
        
        if ~isempty(obj.ste)
            ste = Inf .* ones(obj.volInfo.nvox, k);
            ste(obj.volInfo.wh_inmask, i) = obj.ste(:, i);
        end
        
        if ~isempty(obj.sig)
            sig = zeros(obj.volInfo.nvox, k);
            sig(obj.volInfo.wh_inmask, i) = obj.sig(:, i);
        end
        
    end
    
    
end % function


%  validateattributes Check validity of array.
%     validateattributes(A,CLASSES,ATTRIBUTES) validates that array A belongs
%     to at least one of the specified CLASSES and has all of the specified
%     ATTRIBUTES. If A does not meet the criteria, MATLAB issues a formatted
%     error message.
%  
%     validateattributes(A,CLASSES,ATTRIBUTES,ARGINDEX) includes the
%     position of the input in your function argument list as part of any
%     generated error messages.
%  
%     validateattributes(A,CLASSES,ATTRIBUTES,FUNCNAME) includes the
%     specified function name in generated error identifiers.
%  
%     
%     ATTRIBUTES Cell array that contains descriptions of valid attributes
%                for array A. For example, if ATTRIBUTES = {'real','finite'}
%                A must contain only real and finite values.
%  
%                Supported attributes include:
%     
%                  2d         3d         binary       ndims           nonzero        
%                   <      ncols          nrows      column       nonnegative
%                  <=       size       nonempty        real        decreasing
%                   >        odd      nonsparse      nonnan        increasing
%                  >=      numel         scalar      square     nondecreasing
%                 row       diag        integer      vector     nonincreasing
%                even     finite       positive  scalartext
%                                   
%  
%                Some attributes also require numeric values. For those 
%                attributes, the numeric value or vector must immediately 
%                follow the attribute name string. For example,
%                {'>=', 5, '<=', 10, 'size', [3 4 2]} checks that all
%                values of A are between 5 and 10, and that A is 3-by-4-by-2.
%                