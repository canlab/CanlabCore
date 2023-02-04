function [component_scores eigenmap_obj] = pca(obj, varargin)
% Compute and display principal components for a set of images in obj
%
% :Usage:
% ::
%
%     [component_scores eigenmap_obj] = pca(obj, ['noplot'], ['k', num_components])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Tor Wager, 2023
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
%
% :Inputs:
%
%   **obj:**
%        An image_vector object with multiple images (e.g., an fmri_data
%        object with time series images)
%
% :Optional Inputs:
%   **'noplot':**
%        Suppress plots
%
%   **'k', num_components:**
%        String 'k' followed by number of components to calculate
%        Default = 3, fewer = faster
%
% :Outputs:
%
%   **component_scores:**
%        Score values for each image (rows) on each component (columns)
%
%   **eigenmap_obj:**
%        an image_vector or fmri_data object containing eigenvectors.
%        These correspond to voxel weights on the components.
%        The object contains one image for each component map.
%
% :Examples:
% ::
%
%    obj = load_image_set('emotionreg');
%    [component_scores eigenmap_obj] = pca(obj);
%    % OR: 
%    obj.pca
%
%    See the code of this function for simple ways to create the output
%    plots that are returned by default.
%
%    For time series objects, it may be useful to plot the frequency
%    spectrum of component scores. The code below does this for a run with
%    TR = 1.8 sec.
%    create_figure('fft'); [fftmag, fftfreq, h] = fft_calc(component_scores, 1.8);
%
% :References:
%   There are many publications on PCA, e.g., Herve Abdi's book.
%
% :See also:
%   - pca.m and princomp.m from Matlab
%   - fmri_data.ica, icatb functions.
%

% ..
%    Programmers' notes:
%    Created by Tor Wager, Jan 2023
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

k = 3;
doplot = true;

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'k', 'numcomponents'}, k = varargin{i+1}; varargin{i+1} = [];
            case 'noplot', doplot = false;

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

[eigenmaps, component_scores] = pca(obj.dat', 'Centered', true, 'Economy', true, 'NumComponents', k);

eigenmap_obj = obj;
eigenmap_obj.dat = eigenmaps;

% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
if doplot

    create_figure('components');
    plot_matrix_cols(component_scores);
    set(gca, 'FontSize', 14)
    axis tight
    xlabel('Image number')
    set(gca, 'YColor', 'none')

    for i = 1:k

        comp_obj = get_wh_image(eigenmap_obj, i);
        plot(comp_obj, 'montages');

        canlab_redblue_symmetric_colormap;
        % This function sets the colormap for the current axis to be
        % symmetrical around 0, with positive values shown in red and
        % negative values in blue. It is useful for displaying images and
        % matrices where the sign of the elements and the distribution of 
        % positive and negative values is meaningful. 

        colorbar

        set(gcf, 'Name', sprintf('Component %d', i), 'Tag', sprintf('Component %d', i));

    end

end

end % main function

