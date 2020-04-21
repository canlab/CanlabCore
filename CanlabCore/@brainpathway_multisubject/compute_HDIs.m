% Compute hub disruption indices (HDIs) using BCT functions. See Achard et
% al 2012; Mansour et al 2016.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Yoni Ashar
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
%   **bs:**
%        brainpathway_multisubject object. Will use connectivity.regions.r
%        for computing the HDIs.
%
%   **refgroup:**
%        Two options: 1) a brainpathway_multisubject object to use
%        as a reference group. Will compute BCT measures on refgroup.connectivity.regions.r
%        and take the mean value as a the reference distribution.
%        2) a vector on integers to use for cross-validated estimation of
%        HDIs, where each integer represents one fold. The matlab functions
%        crossvalind and cvpartition can be helpful in generating these
%        indices. The hold-out set is used as the reference distribution
%
% :Optional Inputs:
%   **thresh:**
%        Followed by a link density, ranging from .01 - .99. Will use .1 by default.
%
% :Outputs:
%
%   **bs:**
%        Will save the HDI values in bs.graph_prop, overwriting any
%        existing data in that field
%
%   **refgroup:**
%        With the BCT measures saved in the graph_properties.region field
%
% :Examples:
% ::
%
%    load('brainpathway_data_gsr_censoring_pipeline_nosmoothing.mat');
%    load('paingen_reference_group.mat');
%    
%    [bs, refgroup] = compute_HDIs(bs, bs_paingen)
%
% :References:
%   Achard et al 2012, Mansour et al 2016
%
%
% ..
%    Programmers' notes:
%    Initial commit -- Yoni Ashar, April 2020
% ..
function [bs, refgroup] = compute_HDIs(bs, refgroup, varargin)

thresh = .1; % default 10% link density

allowable_inputs = {'thresh'};

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

%% compute BCT measures on my subjects
print_header(sprintf('Computing BCT measures on %d subjects in test reference group', size(bs.connectivity.regions.r,3)));
bs = bct_toolbox_undirected_graph_metrics(bs, thresh);


%% Compute BCT measures on ref group and regress

% what kind of reference group are we dealing with?

% CV
if isnumeric(refgroup)
    error('CV not implemented yet')
    
% External reference group
elseif isa(refgroup, 'brainpathway_multisubject')
    
    print_header(sprintf('Computing BCT measures on %d subjects in external reference group', size(refgroup.connectivity.regions.r,3)));
    
    %% compute BCT measures on ref group  
    refgroup = bct_toolbox_undirected_graph_metrics(refgroup, thresh);
    
    %% for each subject, regress on mean referece group to get an HDI (i.e., the slope)

    for i=1:size(bs.connectivity.regions.r,3)

        
        % dont compute HDI on global metrics, only nodal
        % TODO: put in subfunction
        mymetrics = bs.graph_properties.regions.Properties.VariableNames;
        metrics = {};
        for m=1:length(mymetrics)
            if numel(bs.graph_properties.regions.(mymetrics{m})(1,:)) > 1
                metrics{end+1} = mymetrics{m};
            end
        end

        % compute the HDIs    
        for m=1:length(metrics) % for each metric

            betas = regress(mean(refgroup.graph_properties.regions.(metrics{m}))', [ones(489, 1) bs.graph_properties.regions.(metrics{m})(i,:)' - mean(refgroup.graph_properties.regions.(metrics{m}))']);

            %figure; scatter(mean(graph_prop_reference_group.(metrics{m})), graph_prop.(metrics{m})(i,:) - mean(graph_prop_reference_group.(metrics{m}))), lsline, title(metrics{m})
            bs.HDIs.regions.(metrics{m})(i) = betas(2);
        end
    end
    
% Unknown
else
    error('Unknown type of reference group')
end


