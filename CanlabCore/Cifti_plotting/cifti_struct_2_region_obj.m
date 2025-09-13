function r = cifti_struct_2_region_obj(cifti_struct, varargin)
% cifti_struct_2_region_obj converts a cifti model structure to a CANlab region object.
% - assumes that grayordinates/voxels are in dimension 1 of the CIFTI
% - returns regions for "models" only with volumetric structures with voxel coordinates, no surfaces
%
% :Optional Inputs:
%
%   **'which_image':** [numeric scalar]
%        Index of the image (map) to render from the cifti file.
%        Default = 1.
%
%   **'verbose'** Print verbose output
%

% Parse variable inputs using inputParser
% --------------------------------------------------------
ARGS = parse_inputs(varargin{:});

fn = fieldnames(ARGS);
for i = 1:length(fn)
    eval([fn{i}, ' = ARGS.(fn{i});']);
end

[model_names, model_types] = listModelNames(cifti_struct);
% see also:
% [surflist, vollist] = cifti_diminfo_dense_get_structures(cifti_struct.diminfo{1});

wh = strcmp(model_types, 'vox'); % only volumetric structures with voxel coordinates, no surfaces

% Check for requirements
filename = which('cifti_struct_dense_extract_surface_data');
if isempty(filename), error('You need Wash U HCP CIFTI tools on your matlab path. See https://github.com/Washington-University/HCPpipelines'); end

% Construct region object
% --------------------------------------------------------

r = [];

for i = 1:length(model_names)

    if wh(i)

        if isempty(r)
            r = construct_region(cifti_struct, i, which_image);

        else

            r(end+1) = construct_region(cifti_struct, i, which_image);

        end

    end

end

end % main function


% ------------------------------------------------------------------------
% Subfunction: parse_inputs
% ------------------------------------------------------------------------
function ARGS = parse_inputs(varargin)
% parse_inputs parses optional input arguments.
%
% :Usage:
% ::
%     ARGS = parse_inputs(optional_name_value_pairs)
%
% :Optional Inputs:
%
%   **'which_image':** [numeric scalar]
%        Index of the image (map) to render from the cifti file.
%        Default = 1.
%
%   **'verbose'** Print verbose output
%
% ------------------------------------------------------------------------
p = inputParser;

addParameter(p, 'which_image', 1, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'verbose', false, @(x) islogical(x) && isscalar(x));

parse(p, varargin{:});
ARGS = p.Results;

end


% ------------------------------------------------------------------------
% Subfunction: listModelNames
% ------------------------------------------------------------------------
function [model_names, model_types] = listModelNames(cifti_struct)
% listModelNames returns a cell array of model names from the cifti structure.

models = cifti_struct.diminfo{1}.models;

numModels = numel(models);
[model_names, model_types] = deal(cell(numModels, 1));

for i = 1:numModels
    model_names{i} = models{i}.struct;
    model_types{i} = models{i}.type;
end

end % listModelNames


% ------------------------------------------------------------------------
% Subfunction: construct_region
% ------------------------------------------------------------------------
function r = construct_region(cifti_struct, indx, which_image)

r = region();

structname = cifti_struct.diminfo{1}.models{indx}.struct;

infostruct = cifti_diminfo_dense_get_volume_structure_info(cifti_struct.diminfo{1}, cifti_struct.diminfo{1}.models{indx}.struct);
% see also
% outinfo = cifti_diminfo_dense_get_volume_all_info(cifti_struct.diminfo{1}, false);

% CIFTI appears to start voxels at 0, so add 1
r.XYZ = cifti_struct.diminfo{1}.models{indx}.voxlist + 1;

r.dim = cifti_struct.diminfo{1}.vol.dims;
r.M = cifti_struct.diminfo{1}.vol.sform;

% r.XYZmm = voxel2mm(r.XYZ, r.M);
r.XYZmm = infostruct.coordlist; % already converted using cifti tools

r.shorttitle = cifti_struct.diminfo{1}.models{indx}.struct;

r.voxSize = abs(diag(r.M(1:3, 1:3)));
r.numVox = size(r.XYZ, 2);
r.center = mean(r.XYZ, 2)';
r.mm_center = mean(r.XYZmm, 2)';
r.descrip1 = cifti_struct.metadata(2).value;
r.title = r.shorttitle;

% r.Z = ones(1, size(r.XYZ, 2)); % placeholder for values that may be added later

% Extract data using HCP tools
datavalues = cifti_struct.cdata(infostruct.ciftilist, which_image);

r.Z = datavalues';

% see also:
% datavalues = cifti_struct_dense_extract_surface_data(cifti_struct, structname, 1);
% datavalues = cifti_struct_dense_extract_volume_structure_data(cifti_struct, structname, true, 1);



% Add data values
% --------------------------------------------------------



end % main construct_region
