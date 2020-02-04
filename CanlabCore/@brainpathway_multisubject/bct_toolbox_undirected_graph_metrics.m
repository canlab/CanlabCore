function bs = bct_toolbox_undirected_graph_metrics(bs)
% bs = bct_toolbox_undirected_graph_metrics(bs)
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


for i = 1:n
    
    [graph_prop, graph_prop_glob] = bct_toolbox_undirected_graph_metrics(r(:, :, i));
    
    vnames = graph_prop.Properties.VariableNames;
    
    bs.graph_properties.regions.metric_names = vnames;
    
    for j = 1:length(vnames)
        
        bs.graph_properties.regions.(vnames{j})(:, i) = graph_prop.(vnames{j});
        
    end
    
    vnames = graph_prop_glob.Properties.VariableNames;
    
    bs.graph_properties.regions.glob_metric_names = vnames;
    
    for j = 1:length(vnames)
        
        bs.graph_properties.regions.(vnames{j})(1, i) = graph_prop_glob.(vnames{j});
        
    end
    
end % subject loop

end % function


