% :Usage:
% ::
%
%     function XYZout = mm2voxel(XYZ, M, [uniqueness options])
%
% :Inputs:
%
%   **M:**
%        VOL must contain M or mat field, which is mat file matrix in SPM
%
%   **XYZ:**
%        can be either 3 rows or 3 columns
%        by transforming to x y z in ROWS, coords in COLUMNS
%        so if you have 3 coordinates, you'd better put the 3 coords
%        in different COLS, with ROWS coding x, y, z!
%
% If a 3rd argument is entered, enter either:
%   1. uniqueness not required; allows repeats
%   2. unique and sorted
%
% Otherwise, voxels are unique and relative order is preserved. 
%
% NB: Unless you require that relative order be preserved, *ALWAYS* set a
%     uniqueness option. They are orders of magnitude faster.

function XYZout = mm2voxel(XYZ, M, varargin)
    XYZout = [];

    if isempty(XYZ)
        return;
    end

    if isstruct(M)
        if isfield(M, 'mat')
            M = M.mat;
        elseif isfield(M, 'M')
            M = M.M;
        else
            error('Unable to identify affine matrix in structure in M\n');
        end
    end


    % transpose XYZ so that x y z are in rows and coords are in cols
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if size(XYZ,1) ~= 3
        XYZ = XYZ';
    end

    XYZout = [XYZ; ones(1, size(XYZ,2))];
    XYZout = (M\XYZout)';
    XYZout(:,4) = [];
    XYZout = round(XYZout);
    XYZout(XYZout == 0) = 1;

    % make sure voxel coords are integers, and there are no repetitions
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % XYZout = unique(round(XYZout),'rows');

    % long method: do not sort (unique.m sorts rows)
    if nargin < 3
        i = 1;
        while i < size(XYZout,1)
            whch = find(ismember(XYZout, XYZout(i,:), 'rows'));
            if length(whch > 1)
                XYZout(whch(2:end),:) = [];
            end
            i = i+1;
        end
    else
        if varargin{1} == 2
            XYZout = unique(round(XYZout),'rows');
        elseif varargin{1} == 1
            % do nothing
        else
            warning('mm2voxel: Enter 1 or 2 for optional argument!')
        end
    end
end
