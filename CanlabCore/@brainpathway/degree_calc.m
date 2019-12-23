
function b2 = degree_calc(b2)
% b2 = degree_calc(b2)
%
% Method for brainpathway
% uses node_clusters to caculate overall, within-cluster(network), and
% between-cluster degree
%
% NOTE: Currently uses region averages (.region), update to use nodes!

node_cluster_numbers = unique(b2.node_clusters);

u = unique(b2.node_clusters);
k = length(u);

[widegree, btdegree, wi_hubs, bt_hubs] = deal(cell(1, k));

for i = 1:k
    
    wh = u(i);  % target network
    
    wh_nodes = b2.node_clusters == wh;
    
    [widegree{i}, btdegree{i}, wi_hubs{i}, bt_hubs{i}] = get_within_between_degree(b2, wh);
    
    
end

b2.graph_properties.regions.within_network_degree = cat(1, widegree{:});
b2.graph_properties.regions.between_network_degree = cat(1, btdegree{:});

b2.graph_properties.regions.within_hubs = cat(1, wi_hubs{:});
b2.graph_properties.regions.between_hubs = cat(1, bt_hubs{:});

% regions x networks matrix, logical
% b2.graph_properties.regions.within_hubs(:, i) = wi_hub_vec;
% b2.graph_properties.regions.between_hubs(:, i) = bt_hub_vec;

b2.graph_properties.regions.descrip = 'weighted degree, mean correlation';

end % main function



function [widegree, btdegree, wi_hubs, bt_hubs] = get_within_between_degree(b2, wh)

wh_nodes = b2.node_clusters == wh;

within_connections = b2.connectivity.regions.r(wh_nodes, wh_nodes);

cross_connections = b2.connectivity.regions.r(wh_nodes, ~wh_nodes);

% ways of cross-connecting:
% 1. simple above threshold
%
% thr = prctile(cross_connections(:), 95);
% [rows, cols] = find(cross_connections > thr);
%
% 2. degree
% connect_degree = sum(cross_connections > thr, 2); % degree of each region, at 95%
% rows = find(connect_degree > prctile(connect_degree, 90));
%
% 3. sum of continuous connections - weighted degree
% note: if we z-score columns we avoid upweighting connections that are
%     themseleves connected; but this may be advantageous.

widegree = mean(within_connections, 2); % degree of each region, at 95%

btdegree = mean(cross_connections, 2); % degree of each region, at 95%

wi_hubs = widegree > prctile(widegree, 90);

bt_hubs = btdegree > prctile(btdegree, 90);


%     Now put these back in overall list
%     -------------------------------------------
% wh_in_target = find(wh_nodes);
% bt_hubs = wh_in_target(bt_hubs);
% wi_hubs = wh_in_target(wi_hubs);

% [bt_hub_vec, wi_hub_vec] = deal(false(wh_nodes, 1));
% 
% bt_hub_vec(bt_hubs) = true;
% wi_hub_vec(wi_hubs) = true;

%     Regions in target network connected to others at >95th percentile
% [bt_hub_vec, wi_hub_vec] = deal(false(size(b2.node_clusters)));
% 
% bt_hub_vec(bt_hubs) = true;
% wi_hub_vec(wi_hubs) = true;

end
