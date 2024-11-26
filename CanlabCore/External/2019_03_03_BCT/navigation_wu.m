function [sr, PL_bin, PL_wei, PL_dis, paths] = navigation_wu(L, D, max_hops)

% Navigation of connectivity length matrix L guided by nodal distance D
%
% % Navigation
% [sr, PL_bin, PL_wei] = navigation_wu(L,D);
% % Binary shortest path length
% sp_PL_bin = distance_bin(L);
% % Weighted shortest path length
% sp_PL_wei = distance_wei_floyd(L);
% % Binary efficiency ratio
% er_bin = mean(mean(sp_PL_bin./PL_bin));
% % Weighted efficiency ratio
% er_wei = mean(mean(sp_PL_wei./PL_wei));
%
% *** Inputs:
%
% L - Weighted/unweighted directed/undirected NxN SC matrix of connection *lengths*
% L(i,j) is the strength-to-length remapping of the connection weight
% between i and j. L(i,j) = 0 denotes the lack of a connection between i
% and j.
%
% D - Symmetric NxN nodal distance matrix (e.g., Euclidean distance between node centroids)
%
% max_hops (optional) - Limits the maximum number of hops of navigation
% paths
%
% *** Outputs:
%
% sr - The success ratio (scalar) is the proportion of node pairs
% successfully reached by navigation.
% 
% PL_bin - NxN matrix of binary navigation path length (i.e., number of hops in
% navigation paths). Infinte values indicate failed navigation paths.
%
% PL_wei - NxN matrix of weighted navigation path length (i.e., sum of connection 
% weights--as defined by C--along navigaiton path). Infinte values indicate failed navigation paths.
%
% PL_dis - NxN matrix of distance-based navigation path length (i.e., sum of connection 
% distances--as defined by D--along navigaiton path). Infinte values indicate failed navigation paths.
%
% paths - NxN cell of nodes comprising navigation paths.
%
% *** Reference: Seguin et al. (2018) PNAS.
% 
% Caio Seguin, University of Melbourne, 2017

    if nargin == 2
        max_hops = length(L);
    end

    N = size(L, 1);
    paths = cell(N);
    PL_bin = zeros(N);
    PL_wei = zeros(N);
    PL_dis = zeros(N);
    
    for i = 1:N
        for j = 1:N
            if (i ~= j)

                curr_node = i;
                last_node = curr_node;
                target = j;
                paths{i,j} = curr_node;
                
                pl_bin = 0;
                pl_wei = 0;
                pl_dis = 0;
                
                while (curr_node ~= target)
                    
                    neighbors = find(L(curr_node,:) ~= 0);
                    
                    [~, min_index] = min(D(target, neighbors));
                    
                    next_node = neighbors(min_index);
                    
                    if isempty(next_node) || next_node == last_node || pl_bin > max_hops

                        pl_bin = Inf;
                        pl_wei = Inf;
                        pl_dis = Inf;
                        break;
                    
                    end
                    
                    paths{i,j} = [paths{i,j} next_node];
                    pl_bin = pl_bin + 1;
                    pl_wei = L(curr_node, next_node) + pl_wei;
                    pl_dis = D(curr_node, next_node) + pl_dis;
                    
                    last_node = curr_node;
                    curr_node = next_node;
                
                end

                PL_bin(i,j) = pl_bin;
                PL_wei(i,j) = pl_wei;
                PL_dis(i,j) = pl_dis;
                
            end
        end
    end
    
    PL_bin(1:N+1:end) = Inf;
    PL_wei(1:N+1:end) = Inf;
    PL_dis(1:N+1:end) = Inf;
    
    sr = 1 - (length(find(PL_bin == Inf)) - N)/(N*N - N);

end
