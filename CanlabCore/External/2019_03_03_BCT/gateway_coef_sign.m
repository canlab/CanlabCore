function [GWpos,GWneg] = gateway_coef_sign(W,Ci,centtype)

%   Gateway coefficient
%
%   [Gpos,Gneg] = gateway_coef_sign(W,Ci,centtype);
%
%   Gateway coefficient is a variant of participation coefficient. Similar
%   to participation coefficient, gateway coefficient measures the
%   diversity of intermodular connections of individual nodes, but this is
%   weighted by how critical these connections are to intermodular
%   connectivity (e.g., if a node is the only connection between it's
%   module and another module, it will have a higher gateway coefficient).
%
%   Inputs:     W,        undirected connection matrix with positive and
%                         negative weights
%
%               Ci,       community affiliation vector
%
%               centtype, centrality measure to use
%                         1 = Node Strength
%                         2 = Betweenness Centrality
%
%   Output:     Gpos,     gateway coefficient for positive weights
%               Gneg,     gateway coefficient for negative weights
%
%   Reference: Vargas ER, Wahl LM. Eur Phys J B (2014) 87:1-10.
%
%   Jeff Spielberg, University of Delaware

%   Modification History:
%   May 2015:  Original (originally adapted from participation_coef_sign.m)
%   July 2018: Bugfix, change in how weighted matrices are handled,
%              improvements for efficiency, additional line documentation

[~,~,Ci]       = unique(Ci);                                        % Remap module indices to consecutive numbers
n              = length(W);                                         % Number of nodes
W(1:(n+1):end) = 0;                                                 % Ensure diagonal is zero
GWpos          = gcoef(W.*(W>0));                                   % Compute gateway coefficient for positive weights
GWneg          = gcoef(-W.*(W<0));                                  % Compute gateway coefficient for negative weights

    function GW = gcoef(W_)
        k    = sum(W_,2);                                           % Compute node strength
        Gc   = (W_~=0)*diag(Ci);                                    % Create neighbor community affiliation matrix
        nmod = max(Ci);                                             % Find # of modules
        ks   = zeros(n,nmod);                                       % Preallocate space
        kjs  = zeros(n,nmod);                                       % Preallocate space
        cs   = zeros(n,nmod);                                       % Preallocate space
        switch centtype                                             % Which centrality measure to use?
            case 1                                                  % Node Strength
                cent = sum(W_,2);
            case 2                                                  % Betweenness Centrality
                L    = weight_conversion(W_,'lengths');
                cent = betweenness_wei(L);
        end
        mcn = 0;                                                    % Set max summed centrality per module to 0
        for i = 1:nmod                                              % For each module
            if sum(cent(Ci==i))>mcn                                 % If current module has a higher sum
                mcn = sum(cent(Ci==i));                             % Reassign value
            end
            ks(:,i) = sum(W_.*(Gc==i),2);                           % Compute the total weight of the connections per node to each module
        end
        for i = 1:nmod                                              % For each module
            if sum(Ci==i)>1                                         % If there is more than 1 node in a module
                kjs(Ci==i,:) = ones(sum(Ci==i),1)*sum(ks(Ci==i,:)); % Compute total module-module connections
                kjs(Ci==i,i) = kjs(Ci==i,i)/2;                      % Account for redundancy due to double counting within-network work weights
            end
        end
        for i = 1:n                                                 % For each node
            if k(i)>0                                               % If node is connected
                for ii = 1:nmod                                     % For each module
                    cs(i,ii) = sum(cent((Ci.*(W_(:,i)>0))==ii));    % Sum of centralities of neighbors of a node within a module
                end
            end
        end
        ksm           = ks./kjs;                                    % Normalize by total connections
        ksm(kjs==0)   = 0;                                          % Account for division by 0
        csm           = cs./mcn;                                    % Normalize by max summed centrality
        gs            = (1-(ksm.*csm)).^2;                          % Calculate total weighting
        GW            = 1-sum((ks.^2)./(k.^2).*gs,2);               % Compute gateway coefficient
        GW(isnan(GW)) = 0;                                          % Account for division by 0
        GW(~GW) = 0;                                                % Set to 0 if no neighbors
    end
end