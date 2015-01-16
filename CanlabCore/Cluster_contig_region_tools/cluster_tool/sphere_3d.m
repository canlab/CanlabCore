% Usage:
% spheres_out = sphere_3d(sphere_centers, r, [unit_length]);
%
% sphere_centers is row or column vector or matrix with rows equal to a
% point in some 3d space, must have 3 columns or rows. If both, it is
% assumed that X Y Z dimensions are specified by different rows, and that
% columns are points. Output is in the same format as the input.
%
% r is a vector of radii with length equal to the number of rows of
% sphere_centers. A scalar input may also be used, resulting in the same
% radii for all spheres.
%
% spheres_out is a matrix with rows equal to all points in the
% spheres defined by sphere_centers and r
%
% unit_length specifies the unit size of the space in each dimension (or a
% scalar for isotropic sizes). default is 1.
%
% sphere_3d assumes that the point 0, 0, 0 exists in the space, and will only
% report points that are integer multiples of unit length distance from
% 0, 0, 0.

function spheres_out = sphere_3d(sphere_centers, r, varargin)
    flip = 0;
    if isempty(sphere_centers) || isempty(r) || (max(size(sphere_centers)) ~= 3 & min(size(sphere_centers))~= 3), error('Error: incorrect inputs'), end
    if size(sphere_centers, 1) ~= 3
        sphere_centers = sphere_centers';
        flip = 1; 
    end

    for i = 1:length(varargin)
        if isnumeric(varargin{i}), 
            if length(varargin{i}) == 1, unit_length = [varargin{i} varargin{i} varargin{i}];
            elseif length(varargin{i}) == 3, unit_length = varargin{i};
            else error('Error: unit_length incorrectly specified')
            end
        else disp(['Warning: do not recognize input ' num2str(i+nargin-length(varargin))])
        end
    end
    if ~exist('unit_length', 'var'), unit_length = [1 1 1]; end

    if length(r) == 1, r(1:size(sphere_centers, 2)) = r; end

    spheres_out = [];

    for i = 1:size(sphere_centers, 2)
        boxmin = space_match(sphere_centers(:,i)-r(i), unit_length, 'force');
        boxmax = space_match(sphere_centers(:,i)+r(i), unit_length, 'force');
        candidates = zeros([3 length(boxmin(1):unit_length(1):boxmax(1))*length(boxmin(2):unit_length(2):boxmax(2))*length(boxmin(3):unit_length(3):boxmax(3))]);
        count = 0;
        for m = boxmin(1):unit_length(1):boxmax(1)
            for n = boxmin(2):unit_length(2):boxmax(2)
                for p = boxmin(3):unit_length(3):boxmax(3)
                    count = count+1;
                    candidates(:,count) = [m n p]';
                end
            end
        end
        for k = 1:size(candidates, 2)
            if sqrt(sum((sphere_centers(:,i)-candidates(:,k)).^2))>r(i)
                candidates(:,k) = NaN;
            end
        end
        [row, col] = find(isnan(candidates));
        col = unique(col);
        if ~isempty(col), candidates(:,col) = []; end
        spheres_out(:,end+1:end+size(candidates, 2)) = candidates;
    end
    spheres_out = space_match(spheres_out, unit_length);
    spheres_out = unique(spheres_out', 'rows')';
    if isempty(spheres_out), spheres_out = space_match(sphere_centers, unit_length, 'force'); end
    if flip, spheres_out = spheres_out'; end
end

