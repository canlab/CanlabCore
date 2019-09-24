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
            
            if isempty(obj.connectivity.nodes)
                create_figure(figtitle, 1, 2);
            else
                create_figure(figtitle, 1, 3);
            end
            
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
          
                subplot(1, 2, 2);

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
            whplot = 2;
            if ~isempty(obj.connectivity.nodes), whplot = 3; end
            
            subplot(1, whplot, whplot);
            title('Ratio of within to between cluster connectivity across Ss');
            
            plot([1 2], [mean(obj.connectivity.regions.within, 2), mean(obj.connectivity.regions.between, 2)])
            hold on
            violinplot({mean(obj.connectivity.regions.within, 2), mean(obj.connectivity.regions.between, 2)}, 'xlabel', {'Within cluster' 'Between cluster'});  
            legend off
            ylabel(char(obj.connectivity_properties.c_fun_han))
         end
       
         % covs should be named fields in subject_metadata. If cell array,
         % first field is of interest, and following fields are covariates
         % of no interest, for computing partial correlations (TODO: not
         % implemented yet)
         function plot_qcfc(obj, covs)
             
             % for convenience
            r = obj.connectivity.regions.r;

            % get indices for the lower triangle, excluding diagonal
            [row, col] = find(triu(r(:,:,1), 1));

            % for each edge, 
            for i=1:length(row)
                % compute compute QC FC
                qcfc(i) = corr(squeeze(r(row(i),col(i),:)), obj.subject_metadata.(covs));
            end

            figure; histogram(qcfc)
            title(covs)
            
         end          	
    end
    
end