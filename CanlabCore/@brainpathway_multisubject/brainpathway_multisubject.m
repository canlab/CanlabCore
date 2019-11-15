% This object class is for multisubject (group level) brainpathway analyses
% It is a child of (inherits from) brainpathway; see documentation there
% for more info. The general idea is that fields in this object will have
% an additional dimension for subject. For example, in brainpathway,
% connectivity is k x k. Here, connectivity is k x k x n. 
%
% The subject_metadata table includes subject ID (string), subject ID
% (integer), and any other meta data of interest (condition, avg head
% motion, age, etc.)
%
% For now, drop voxel-level data when working at the group level
%
% TODO: expand documentation

% if inherit, node_dat and region_dat must be single arrays, which is
% possible to honor in the multi_subject object but perhaps not ideal. For
% now, inherit and don't deal w/ node_dat and region_dat
classdef brainpathway_multisubject < brainpathway
   
    properties
        
        % for now, don't define voxel data here
        % skip for now b/c inheriting
        %node_dat  (1, :) cell; % one cell per subj
        %region_dat (1, :) cell; % one cell per subj
        
        subject_metadata table;
        
    end
    
    
    methods

        function obj = add_subject(obj, new_subj)

            validateattributes([new_subj], {'brainpathway'}, {'nonempty'});
            
            % bring along node and region data; must deal w/ differing num
            % of timepoints for each S, so save each S in a cell array 
            % SKIP FOR NOW B/C CAN'T BE CELL ARRAY
            %obj.node_dat{end+1} = new_subj.node_dat;
            %obj.region_dat{end+1} = new_subj.region_dat;
            
            % if first subj added, grab it all
            if isempty(obj.connectivity.regions)
                obj.connectivity.regions = new_subj.connectivity.regions;
            else % for later Ss, append data
                obj.connectivity.regions.r(:,:,end+1) = new_subj.connectivity.regions.r;
                obj.connectivity.regions.p(:,:,end+1) = new_subj.connectivity.regions.p;
                
                if isfield(new_subj.connectivity.regions, 'within')
                    obj.connectivity.regions.within(end+1, :) = new_subj.connectivity.regions.within;
                    obj.connectivity.regions.between(end+1, :) = new_subj.connectivity.regions.between;
                    obj.connectivity.regions.avg_between_over_within(end+1) = new_subj.connectivity.regions.avg_between_over_within;
                end
            end
            
            row = height(obj.subject_metadata) + 1;
            obj.subject_metadata.id(row) = row;
            
            % median corr in the lower diagonal of corr matrix
            obj.data_quality.median_corr(row) = median(nonzeros(tril(new_subj.connectivity.regions.r)));
            
        end
    
         function plot_connectivity(obj, varargin)
        % Takes any optional input arguments to plot_correlation_matrix
        
            input_args = varargin;
            
            % take mean in 3rd dimension to avg over Ss
            S = struct('r', mean(obj.connectivity.regions.r, 3), 'p', mean(obj.connectivity.regions.p, 3), 'sig', mean(obj.connectivity.regions.p < 0.05, 3));
            
            Xlabels = format_strings_for_legend(obj.region_atlas.labels);
                    
            figtitle = 'brainpathway_connectivity_view_multisubject';
            
            create_figure(figtitle, 2, 2);
            
            k = size(S.r, 1);
            % Circle-plot display and text are automatically suppressed for
            % k > 50 and 15, respectively, in plot_correlation_matrix
            
            if k > 50
                % Plot without text labels
                OUT = plot_correlation_matrix(S, 'nofigure', varargin{:});
                
            else
                
                OUT = plot_correlation_matrix(S, 'nofigure', ...
                    'var_names', Xlabels, varargin{:});
                
            end
            
            num_nodes = length(obj.region_indx_for_nodes);
            if num_nodes == num_regions(obj.region_atlas)
                
                node_labels = Xlabels;
                
            else node_labels = {};
            end
            
            title('Region connectivity')
            
            if ~isempty(obj.connectivity.nodes)
          
                subplot(2, 2, 2);

                % mean in 3rd dimension to avg over Ss
                S = struct('r', mean(obj.connectivity.nodes.r, 3), 'p', mean(obj.connectivity.nodes.p, 3), 'sig', mean(obj.connectivity.nodes.p < 0.05, 3));

                if k > 50
                    % Plot without text labels
                    OUT = plot_correlation_matrix(S, 'nofigure', varargin{:});

                else

                    OUT = plot_correlation_matrix(S, 'nofigure', ...
                        'var_names', node_labels, varargin{:});

                end

                title('Node connectivity')

            end
            
            % add histogram of avg within/between
            if isfield(obj.connectivity.regions, 'within')
                whplot = 4;
                if ~isempty(obj.connectivity.nodes), whplot = 3; end

                subplot(2, 2, whplot);
                title('Ratio of within to between cluster connectivity across Ss');

                plot([1 2], [mean(obj.connectivity.regions.within, 2), mean(obj.connectivity.regions.between, 2)])
                hold on
                violinplot({mean(obj.connectivity.regions.within, 2), mean(obj.connectivity.regions.between, 2)}, 'xlabel', {'Within cluster' 'Between cluster'});  
                legend off
                ylabel(char(obj.connectivity_properties.c_fun_han))
            end
            
            % plot median connectivity value for each S
            subplot(2,2,3)
            histogram(median(obj.connectivity.regions.r, [1 2]));
            title('median connectivity values across Ss')
         end
       
         % Computes pearson correlation for each edge between connectivity
         % and given covariates, across subjects. Covariates can include
         % subject-level measures of motion (e.g., mean FD), in which case
         % correlations between motion and edge strength indicate that
         % motion artifacts are present in the data -- higher motion
         % subjects typically exhibit stronger connectivity (e.g., see
         % Parkes et al 2018 for one example).
         %
         % Input:  covs is a char array or cell array containing fields in subject_metadata. If cell array,
         % first field is of interest, and following fields are covariates
         % of no interest, and partial correlations are computed.
         function qcfc = plot_qcfc(obj, covs)
             
             % for convenience
            r = obj.connectivity.regions.r;

            % get indices for the lower triangle, excluding diagonal
            [row, col] = find(triu(r(:,:,1), 1));
            fprintf('Computing %d correlations, this might take a minute...\n', length(row))
            % for each edge, 
            for i=1:length(row)
                % compute compute QC FC                
                if iscell(covs)
                    qcfc(i) = partialcorr(squeeze(r(row(i),col(i),:)), obj.subject_metadata.(covs{1}), obj.subject_metadata{:,covs{2:end}}, 'Rows', 'pairwise');
                else
                    qcfc(i) = corr(squeeze(r(row(i),col(i),:)), obj.subject_metadata.(covs), 'Rows', 'pairwise');
                end
            end

            create_figure('qcfc'); histogram(qcfc)
            
            if iscell(covs)
                title(sprintf('Partial corrs. between edge connectivity and %s, across Ss', covs{1}))
            else
                title(sprintf('Corrs. between edge connectivity and %s, across Ss', covs))
            end
            
            fprintf('Mean: %3.2f\tMedian: %3.2f\tSD: %3.2f\n', mean(qcfc), median(qcfc), std(qcfc));
         end          	
    end
    
end