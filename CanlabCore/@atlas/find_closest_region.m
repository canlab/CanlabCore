function out = find_closest_region(atl, input_xyz_coord)
% Find the closest atlas region to an input [x, y, z] mm coordinate
%
% First line: One-line summary description of function
%
% :Usage:
% ::
%
%     out = find_closest_region(atlas_obj, input_xyz_coord)
%
%
% :Inputs:
%
%   **atlas_obj:**
%        An atlas-class object
%
%   **input_xyz_coord:**
%        A row vector with x, y, z coordinates in mm
%
% :Outputs:
%
%   **out:**
%        A structure containing:
%
%
% Example:
% Find the Glasser 2016 parcel containing a coordinate and display it
% --------------------------------------------------------------------
% atl = load_atlas('glasser');
% input_xyz_coord = [-18 12 -22];
% out = find_closest_region(atl, input_xyz_coord);
% disp(out)
% orthviews(out.closest_region_obj);
% spm_orthviews('reposition', input_xyz_coord)
%
% Note: requires dist.m from Matlab neural net toolbox to work correctly.
% It is very fast.
% There are many replacements that could be added, but this is the current
% implementation.

% ..
%     Author and copyright information:
%
%     Copyright (C) 2021 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..


out = struct();

if ~isrow(input_xyz_coord)
    input_xyz_coord = input_xyz_coord';
end

if ~isrow(input_xyz_coord)
    error('Input a row vector with [x y z] coordinates in mm');
end

% transform all coordinates in object into mm
mmcoords = voxel2mm(atl.volInfo.xyzlist', atl.volInfo.mat);

% find closest
out.d = dist(input_xyz_coord , mmcoords);

[out.mind, out.wh] = min(out.d);

out.region_number = atl.dat(out.wh);

if length(unique(out.region_number)) > 1
    warning('More than one region is equally close to input coordinate.');
end


out.closest_label = atl.labels{out.region_number};

out.closest_region_obj = select_atlas_subset(atl, out.region_number);

end % function

