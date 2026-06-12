function [dat, tbl] = assign_vals(atl, varargin)
% assign_vals Assign numeric values to atlas regions and create an fmri_data object.
%
% Take an atlas object and a vector of values associated with named
% regions, and produce an fmri_data object whose voxel values are the
% assigned region value, plus a table mapping each region name to its
% assigned value. Region names are matched against the atlas region
% shorttitles (as produced by atlas2region).
%
% :Usage:
% ::
%
%     [dat, tbl] = assign_vals(atl, 'reg_names', reg_names, 'vals', vals, 'sort', true)
%
% :Inputs:
%
%   **atl:**
%        An atlas-class object (with fields .labels, .dat, etc.).
%
% :Optional Inputs:
%
%   **'reg_names':**
%        Cell array (or string array) of region names to assign values to.
%        Names must match atlas region shorttitles. Default: atl.labels.
%
%   **'vals':**
%        Numeric vector of values to assign, the same length as reg_names.
%        Default: zeros(numel(atl.labels), 1).
%
%   **'sort':**
%        Logical flag to sort output table by value in descending order.
%        Default: true.
%
% :Outputs:
%
%   **dat:**
%        fmri_data object with each region's voxel values set to the
%        assigned value.
%
%   **tbl:**
%        MATLAB table of region shorttitles and assigned values.
%
% :Examples:
% ::
%
%     [dat, tbl] = assign_vals(atl, 'reg_names', {'Reg1' 'Reg2'}, ...
%                              'vals', [1 2], 'sort', true);
%
% :See also:
%   - atlas2region
%   - region2fmri_data
%
% ..
%    Author: Michael Sun, Ph.D. 4/23/2025
% ..

    % Parse optional inputs
    parser = inputParser;
    parser.KeepUnmatched = true;

    addParameter(parser, 'reg_names', atl.labels, @(x) iscell(x) || isstring(x));
    addParameter(parser, 'vals', repmat(0, numel(atl.labels), 1), @isnumeric);
    addParameter(parser, 'sort', true, @islogical);
    
    parse(parser, varargin{:});
    reg_names = cellstr(parser.Results.reg_names);  % ensure cell array of char
    vals = parser.Results.vals;
    do_sort = parser.Results.sort;


    % Validate input lengths
    if numel(reg_names) ~= numel(vals)
        error('reg_names and vals must be the same length.');
    end


    % ------------------------------
    % Convert atlas to region object
    % ------------------------------
    atl_r = atlas2region(atl);
    counts = zeros(numel(atl_r), 1);

    % ------------------------------
    % Assign values to regions
    % ------------------------------
    for i = 1:numel(atl_r)
        idx = find(strcmp(reg_names, atl_r(i).shorttitle));
        if ~isempty(idx)
            atl_r(i).Z = repmat(vals(idx), 1, size(atl_r(i).Z, 2));
            counts(i) = vals(idx);
        else
            atl_r(i).Z = zeros(1, size(atl_r(i).Z, 2));
            counts(i) = 0;
        end
    end

    % ------------------------------
    % Create output table
    % ------------------------------
    region_labels = {atl_r.shorttitle}';
    tbl = table(region_labels, counts, 'VariableNames', {'Region', 'Value'});

    if do_sort
        tbl = sortrows(tbl, 'Value', 'descend');
    end

    % ------------------------------
    % Convert region object back to fmri_data
    % ------------------------------
    dat = region2fmri_data(atl_r, fmri_data(atl));

end