function bs = bct_toolbox_undirected_graph_metrics(bs, thresh, varargin)
% bs = bct_toolbox_undirected_graph_metrics(bs)
%
% thresh = proportional link density threshold; 0 to 1 value. (.1 is a common value)
%
% Method for brainpathway_multisubject that extracts graph metrics using
% the Sporns lab (et al.) BCT toolbox.
% Uses bct_toolbox_undirected_graph_metrics
% 
% Stores results in bs.graph_properties.regions
% - Nodes x subjects matrix for each graph property
% - Can do group stats/t-tests/etc on (matrices)'
%
% Examples:
% ------------------------------------------------------------------
% Load brainpathway_data_gsr_censoring_pipeline.mat % from OLP/Yoni
% bs = bct_toolbox_undirected_graph_metrics(bs);
% mean_degree = nanmean(bs.graph_properties.regions.degree')'; % mean degree for each node
%
% % Plot montage of mean degree on brain
% fmri_dat_obj = brainpathway2fmri_data(bs, bs.graph_properties.regions.degree);
% montage(mean(fmri_dat_obj));


r = double(bs.connectivity.regions.r);

n = size(r, 3); % subjects
n_nodes = size(r, 1); % nodes

bs.graph_properties.regions = table(); % clear out -- overwrite

fprintf('Completed  subject        ');

warning off
for i = 1:n
    
    [graph_prop, graph_prop_glob] = bct_toolbox_undirected_graph_metrics(r(:, :, i), thresh, varargin{:});
    
    % nodal properties
    vnames = graph_prop.Properties.VariableNames;    
    for j = 1:length(vnames)
        
        bs.graph_properties.regions.(vnames{j})(i, 1:n_nodes) = graph_prop.(vnames{j});
        
    end
    
    % global properties
    vnames = graph_prop_glob.Properties.VariableNames;    
    for j = 1:length(vnames)
        
        bs.graph_properties.regions.(vnames{j})(i, 1) = graph_prop_glob.(vnames{j});
        
    end
    
    % print to update status
    fprintf('\b\b\b\b\b\b%5d\n',i) 

end % subject loop
warning on

end % function


    


