function [stats, ntmaps, hh, hhfill, table_group, multcomp_group] = hansen_neurotransmitter_maps(fmri_data_obj, varargin)
% Apply a set of PET neurotransmitter-related maps from Hansen 2022 Nature Neuroscience
% - Uses cosine similarity metric on gray matter masked images
%
% :Usage:
% ::
%
%     [stats, ntmaps, hh, hhfill, table_group, multcomp_group] = hansen_neurotransmitter_maps(fmri_data_obj, [optional inputs])
%
% :Inputs:
%
%   **fmri_data_obj:**
%        An fmri_data object with one or more images; can be a signature
%
% :Optional Inputs:
%   **colors:**
%        Followed by plot colors in cell array
%
%   **nofigure:**
%        Do not create a new figure; use existing
%
%   **noplot:**
%        Do not plot polar plot, return tables only
%
%   **correlation:**
%        Use correlation similarity metric instead of cosine similarity
%
%   **dofixrange:**
%        Followed by range of values [min max] for plot (see
%        image_similarity_plot)
%   
%    **doAverage:**
%      Plot average correlation from the input images
% :Outputs:
%
%   **stats:**
%        Stats structure with associations
%
%   **ntmaps:**
%        fmri_data object with neurotransmitter maps and meta-data, in
%        sorted order
%
%   **'hh'**
%        line handles
%
%   **'hhfill'**
%        fill handles
%
%    **table_group**
%              multiple one-way ANOVA tables (one for each
%              spatial basis) with group as column factor (requires
%              'average' to be specified)
%  
%    **multcomp_group**
%              mutiple comparisons of means across groups, one output
%              cell for each spatial basis, critical value determined
%              by Tukey-Kramer method (see multcompare)
%
% :Examples:
% ::
%
% % Load and plot profile for a negative affect neuromarker from Chang et al. 2015, Plos Biology:
% % ------------------------------------------------------------------
% pines = load_image_set('pines');
% [stats, ntmaps] = hansen_neurotransmitter_maps(pines);
%
% figure;
% [stats, ntmaps] = hansen_neurotransmitter_maps(pines, 'nofigure', 'colors', [0 0 1]);
% [stats, ntmaps] = hansen_neurotransmitter_maps(pines, 'noplot');
%
% % Load a drug craving neuromarker from Koban et al. 2022, Nature
% % Neuroscience, and plot the profiles for drug-only and drug+food signatures on top of
% % one another:
% % ------------------------------------------------------------------
% ncs_all = load_image_set('ncs'); % craving
% ncsdrug = get_wh_image(ncs_all, 2);
% ncs = get_wh_image(ncs_all, 1);
% [stats, ntmaps, hh, fillh, table_group] = hansen_neurotransmitter_maps(ncsdrug);
% [stats, ntmaps, hh, fillh, table_group] = hansen_neurotransmitter_maps(ncs, 'colors', {[0 0 1]}, 'nofigure');

% :References:
%   Hansen et al. 2022, Nature Neuroscience
%
% :See also:
%   - image_similarity_plot
%

%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2023 Mijin Kwon, Tor Wager and Ke Bo
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

% ..
%    Programmers' notes:
%    Tor: Created from Mijin's script -- plus resort rows of image names and
%    fullpath fields

%    Ke:Add function to input group of image_obj and plot mean correlation.
%    e.g group of bootsrap sample. In this case, the error bar is displayed
%    by standard deviation.
% ..

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

colors = [1 0 0];
dofigure = true;
doplot = true;
similarity_metric = 'corr';
dofixrange = [];

doAverage=0;
% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'colors' 'doplot' 'similarity_metric' 'dofixrange'};

keyword_inputs = {'noplot' 'nofigure' 'cosine_similarity','doAverage'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs

                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);

            case keyword_inputs
                % Skip deal with this below
            
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% 2nd pass: Keyword inputs. These supersede earlier inputs
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'noplot'
                doplot = false;

            case 'nofigure'
                dofigure = false;

            case 'cosine_similarity'
                similarity_metric = 'cosine_similarity';
            case 'doAverage'
                doAverage=1;

        end
    end
end

% -------------------------------------------------------------------------
% MAIN FUNCTION
% -------------------------------------------------------------------------

ntmaps = load_image_set('hansen22');

ntmaps = reorder_and_add_metadata(ntmaps);

% These are already gray-matter masked in repo, but make sure:
ntmaps = apply_mask(ntmaps, which('gray_matter_mask.nii'));

% This may not be masked...so mask with gray matter:
fmri_data_obj = apply_mask(fmri_data_obj, which('gray_matter_mask.nii'));

if dofigure
    create_figure('Neurotransmitter polar plot')
end

if doplot
    if ~iscell(colors), colors = {colors}; end

    if doAverage==1
        if isempty(dofixrange)
             [stats, hh, hhfill, table_group, multcomp_group] = image_similarity_plot(fmri_data_obj, 'mapset', ntmaps, similarity_metric, 'plotstyle', 'polar', 'networknames', ntmaps.metadata_table.target, 'colors', colors, 'nofigure', 'average','Error_STD');
        else

              [stats, hh, hhfill, table_group, multcomp_group] = image_similarity_plot(fmri_data_obj, 'mapset', ntmaps, similarity_metric, 'plotstyle', 'polar', 'networknames', ntmaps.metadata_table.target, 'colors', colors, 'nofigure', 'dofixrange', dofixrange,'average','Error_STD');
        end
    else
        if isempty(dofixrange)
    
            [stats, hh, hhfill, table_group, multcomp_group] = image_similarity_plot(fmri_data_obj, 'mapset', ntmaps, similarity_metric, 'plotstyle', 'polar', 'networknames', ntmaps.metadata_table.target, 'colors', colors, 'nofigure');
    
        else % we have fixed range
    
            [stats, hh, hhfill, table_group, multcomp_group] = image_similarity_plot(fmri_data_obj, 'mapset', ntmaps, similarity_metric, 'plotstyle', 'polar', 'networknames', ntmaps.metadata_table.target, 'colors', colors, 'nofigure', 'dofixrange', dofixrange);
    
        end
    end

else
    [stats, hh, hhfill, table_group, multcomp_group] = image_similarity_plot(fmri_data_obj, 'mapset', ntmaps, similarity_metric, 'noplot');

end

end % main function



% ---------------------------------------------------------------
% Subfunctions
% ---------------------------------------------------------------

function ntmaps = reorder_and_add_metadata(ntmaps)

nrows = height(ntmaps.metadata_table);

transmitter = cell(nrows, 1);
ntmaps.metadata_table = addvars(ntmaps.metadata_table, transmitter,'After','target');
ntmaps.metadata_table.Properties.VariableNames{2} = 'transmitter';

for i = 1:nrows

    if ismember(ntmaps.metadata_table.target{i}, {'5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT'})
        ntmaps.metadata_table.transmitter{i} = 'serotonin';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'CB1'})
        ntmaps.metadata_table.transmitter{i} = 'cannabinoid';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'D1', 'D2', 'DAT'})
        ntmaps.metadata_table.transmitter{i} = 'dopamine';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'GABAa', 'GABAabz'})
        ntmaps.metadata_table.transmitter{i} = 'GABA';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'H3'})
        ntmaps.metadata_table.transmitter{i} = 'histamine';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'M1', 'a4b2', 'VAChT'})
        ntmaps.metadata_table.transmitter{i} = 'acetylcholine';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'MOR'})
        ntmaps.metadata_table.transmitter{i} = 'opioid';
    end
    if ismember(ntmaps.metadata_table.target{i}, {'NET'})
        ntmaps.metadata_table.transmitter{i} = 'norepinephrine';
    end

    if ismember(ntmaps.metadata_table.target{i}, {'mGluR5'})
        ntmaps.metadata_table.transmitter{i} = 'glutamate';
    end
end

% display(ntmaps.metadata_table)

[ntmaps.metadata_table, wh_indx_all] = sortrows(ntmaps.metadata_table,'transmitter','ascend');

ntmaps.dat = ntmaps.dat(:, wh_indx_all);

% Save original image names in sorted order
ntmaps.metadata_table = addvars(ntmaps.metadata_table, ntmaps.image_names(wh_indx_all, :), 'NewVariableName', {'image_names'}, 'Before', 'modeling_notes');

% re-sort
ntmaps.fullpath = ntmaps.fullpath(wh_indx_all, :);

end % end subfunction



