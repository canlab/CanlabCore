%% load data

load('FisherIris.mat')
rng('default');

%% Run umap

umap_coords = run_umap(meas);

% see also 'graphic' cluster output

%%
f1 = create_figure('umap on iris dataset', 1, 3); 

plot(umap_coords(:, 1), umap_coords(:, 2), 'ko');

xlabel('umap(1)'); ylabel('umap(2)');
title('UMAP with colors indicating true classes'); 

[indic, names, condf] = string2indicator(species);
colors = scn_standard_colors(length(names));

for i = 1:length(names)
    
    wh = condf == i;
    
    plot(umap_coords(wh, 1), umap_coords(wh, 2), 'ko', 'MarkerFaceColor', colors{i});
    
end

drawnow

%% Cluster data in derived umap space and re-plot

% note: clusterIdentifiers from run_umap didn't produce great results for me (tor)

% note: you really want enough (> 2) umap dimensions for this to be meaningful.
% Run again with 3 umap dims
% Will run with same number of dims as original dimensions, but not sure
% yet if this is a good idea.

[umap_coords3d, umap, umap_clusterIdentifiers, extras] = run_umap(meas, 'cluster_output', 'numeric', 'n_components', 3);

clusterIdentifiers = clusterdata(umap_coords3d, 'linkage','ward','savememory','on','maxclust', 4);

%% re-plot

figure(f1)
subplot(1, 3, 2)

plot(umap_coords(:, 1), umap_coords(:, 2), 'ko');

xlabel('umap(1)'); ylabel('umap(2)');
title('UMAP with colors indicating estimated clusters'); 

n = length(unique(clusterIdentifiers));

colors = scn_standard_colors(n);

for i = 1:n
    
    wh = clusterIdentifiers == i;
    
    plot(umap_coords(wh, 1), umap_coords(wh, 2), 'ko', 'MarkerFaceColor', colors{i});
    
end

drawnow

%%
figure(f1)
subplot(1, 3, 3)

plot(umap_coords3d(:, 1), umap_coords3d(:, 2), 'ko');

xlabel('umap(1)of 3'); ylabel('umap(2) of 3');
title('UMAP with colors indicating estimated clusters'); 

n = length(unique(clusterIdentifiers));

colors = scn_standard_colors(n);

for i = 1:n
    
    wh = clusterIdentifiers == i;
    
    plot(umap_coords3d(wh, 1), umap_coords3d(wh, 2), 'ko', 'MarkerFaceColor', colors{i});
    
end

drawnow
