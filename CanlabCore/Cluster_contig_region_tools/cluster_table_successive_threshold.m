function cl = cluster_table_successive_threshold(cl,varargin)
% Cluster table of a cell array of clusters cl{1} cl{2} etc.
% Prints table of cl{1} and then any additional regions in cl{2:n} that are
% not within 10 mm of a previously printed cluster
%
% Also: merges clusters in set within 10 mm
% 
% Table titles are hard-coded to be consistent with meta-analysis toolbox right now
%
% :Usage:
% ::
%
%    cl = cluster_table_successive_threshold(cl,[sizethr])
%
% :Example:
% ::
%
%    % Print a series of tables with custom fields:
%    cl = cluster_table_successive_threshold(cl,5,'myfield1','myfield2')
%
% Use *merge_nearby_clusters* and *subclusters_from_local_max*
% or some other way to get clusters appropriately separated and distanced
% before running.
% see Meta_cluster_tools for code to run this for meta-analysis.
%
% ..
%    tor wager, oct/nov 06
% ..

higher_thr = [];
mind = 10;          % distance
showsubpeaks = 0; 	% applies to 2nd-nth cells of clusters

if nargin < 2, sizethr = 0; else sizethr = varargin{1}; end

doextrafields = 0;
if nargin > 2
    doextrafields = 1;
    % we have additional fields in clusters to print in table
    myfields = varargin(2:end);
end

    

% merge nearby clusters
% % for i = 1:length(cl)
% %     if ~isempty(cl{i})
% %         cl{i} = merge_nearby_clusters(cl{i},mind);
% %     end
% % end
% % 
% % if doextrafields
% %     sum_merged_vals;
% % end
    
fprintf(1,'Height\n');
if doextrafields
    base = 'cluster_table(cl{1},showsubpeaks,0';
    estr = build_call(base,myfields);
    eval(estr);
else
    cluster_table(cl{1},1,0);
end

fprintf(1,'\n');


if isempty(cl{1})
    % do nothing
else
% %     subc = subclusters_from_local_max(cl{1}, mind);
% %     higher_thr = subc;
    higher_thr = cl{1};
end

strs = {'Stringent' 'Medium' 'Lenient'};

for i = 2:length(cl)
    
    fprintf(1,'Additional regions at Extent: %s and size >= %3.0f\n',strs{i-1},sizethr);
    % %     subc = subclusters_from_local_max(cl{i}, 10);

    if ~isempty(cl{i})
        subc = cl{i};

        if isempty(higher_thr)
            dosubctable = showsubpeaks;    % subclusters in table, this is first table...
            outside_range = ones(1,length(subc));
        else
            dosubctable = 0;
            [close_enough,outside_range,nearest_distance,closest_cluster] = cluster_close_enough(higher_thr,subc,mind);
        end

        % size threshold
        sz = cat(1,subc.numVox); sz = sz' >= sizethr;

        if doextrafields
            base = 'cluster_table(subc(outside_range & sz),dosubctable,0';
            estr = build_call(base,myfields);
            eval(estr);
        else
            cluster_table(subc(outside_range & sz),dosubctable,0);
        end

        fprintf(1,'\n');

        higher_thr = merge_clusters(higher_thr,subc);

        % save for output, to plot and stuff
        cl{i} = subc(outside_range & sz);
    end
end


% nested functions
% %     function sum_merged_vals
% %         disp('Summing results of extra fields for merged clusters.')
% %         disp('Appropriate for meta-analysis, but maybe not for other applications.')
% %         for i = 1:length(cl)
% %             for j = 1:length(cl{i})
% %                 for k = 1:length(myfields)
% %                     cl{i}(j).(myfields{k}) = sum(cl{i}(j).(myfields{k}));
% %                 end
% %             end
% %         end
% %     end
    

end  % main function




function estr = build_call(estr,myfields)
    % build table function call
%estr = 'cluster_table(subc(outside_range & sz),dosubctable,0';

for i = 1:length(myfields)
    estr = [estr ',''' myfields{i} ''''];
end
estr = [estr ');'];

end


    
