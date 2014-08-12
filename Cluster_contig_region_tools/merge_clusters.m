% subclusters = merge_clusters(clusters_to_match, subclusters_to_change)
% MERGE_CLUSTERS - function for synchronizing the field list of cluster structures
%   and merging them
function outcl = merge_clusters(c2m, subcl)
    if isempty(c2m)
        disp('Merge clusters: Clusters to match is empty.');
        outcl = subcl;
        return
    end

    if isempty(subcl)
        outcl = c2m;
        return
    end

    if ~iscell(subcl)
        s{1} = subcl;
        subcl = s;
    end
    outcl = c2m;

    for i = 1:length(subcl)
        for j = 1:length(subcl{i})
            % make all field names match.
            newcl{i}(j) = c2m(1);
            N = fieldnames(newcl{i}(j));
            for NN = 1:length(N)
                eval(['newcl{i}(j).' N{NN} ' = [];'])
            end

            %N = fieldnames(clusters(i));
            for NN = 1:length(N)
                str = ['newcl{i}(j).' N{NN} ' = subcl{i}(j).' N{NN} ';'];
                try
                    eval(str)
                catch
                    disp(['Making ' N{NN}])
                    str = ['newcl{i}(j).' N{NN} ' = [];'];
                    eval(str)
                end
            end
        end

        outcl = ([outcl, newcl{i}]);
    end
end
