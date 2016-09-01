function [sim_matrix, layer2names, layer2colors] = riverplot_reorder_matrix(sim_matrix, layer2names, layer2colors)


% sim_matrix
% names
% colors

% pca, sort scores ascending based on first component
% [coeff, score, latent, tsquare] = pca(sim_matrix, 'NumComponents', 1);

c = clusterdata(sim_matrix,'linkage','ward','savememory','on', 'maxclust', round(size(sim_matrix, 1) ./2));
[cs, indx] = sort(c, 1, 'ascend');

sim_matrix = sim_matrix(indx, :);

%layer2 = layer2(indx);
layer2names = layer2names(indx);
layer2colors = layer2colors(indx);

end
