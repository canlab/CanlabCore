function [dat, tbl] = assign_vals(atl, varargin)
% ASSIGN_VALS Assigns numeric values to atlas regions and creates an fmri_data object.
%
%   [dat, tbl] = assignvals2atlas(atl, 'reg_names', reg_names, 'vals', vals, 'sort', true)
%
%   Inputs:
%     atl        - An atlas structure (with fields like .labels, .data)
%     reg_names  - Cell array of region names (must match atlas region shorttitles)
%     vals       - Numeric values (same length as reg_names)
%     sort       - (Optional) Logical flag to sort output table by value (default = true)
%
%   Outputs:
%     dat        - fmri_data object with assigned region values
%     tbl        - Table of regions and assigned values
%
% Author: Michael Sun, Ph.D. 4/23/2025

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