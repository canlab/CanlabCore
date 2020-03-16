% Select a subset of subjects
% wh_keep: logical array size, size n_subjects x 1
function b2 = get_wh_subjects(bs, wh_keep)
  
    if ~islogical(wh_keep), error('wh_keep must be a logical array'), end

    fprintf('Creating a copy of brainpathway_multisubject object\n');

    % Need to create a new object first with the same atlas
    b2 = brainpathway_multisubject(bs.region_atlas);          % Construct a brainpathway object from an atlas object

    myfields = fieldnames(bs);  %{'voxel_dat' 'node_dat' 'region_dat' 'network_dat' 'partition_dat'
    myfields(strcmp(myfields, 'region_atlas')) = [];

    for i = 1:length(myfields)

        b2.(myfields{i}) = bs.(myfields{i});
        
        % recursive copy
        if isstruct(bs.(myfields{i}))
            mysubfields = fieldnames(bs.(myfields{i}));
            for j = 1:length(mysubfields)
                b2.(myfields{i}).(mysubfields{j}) = bs.(myfields{i}).(mysubfields{j});
            end
            
            if ~isempty(bs.(myfields{i})) & isstruct(bs.(myfields{i}).(mysubfields{j}))
                mysubsubfields = fieldnames(bs.(myfields{i}).(mysubfields{j}));
                for k = 1:length(mysubsubfields)
                    b2.(myfields{i}).(mysubfields{j}).(mysubsubfields{k}) = bs.(myfields{i}).(mysubfields{j}).(mysubsubfields{k});
                end
            end
        end
    end

    % select subjects
    b2.subject_metadata(~wh_keep,:) = [];
    b2.connectivity.regions.r(~wh_keep,:) = [];
    b2.connectivity.regions.p(~wh_keep,:) = [];
    
    if ~isempty(b2.connectivity.nodes)
        b2.connectivity.nodes.r(~wh_keep,:) = [];
        b2.connectivity.nodes.p(~wh_keep,:) = [];
    end
end