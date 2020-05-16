% Select a subset of subjects
% wh_keep: logical array size, size n_subjects x 1
% WARNING: ONLY PARTIALLY IMPLEMENTED AT THE MOMENT. SEE CODE FOR DETAILS
%
% Yoni Ashar, March 2020
function bs2 = get_wh_subjects(bs, wh_keep)
  
    if ~islogical(wh_keep), error('wh_keep must be a logical array'), end

    fprintf('Creating a copy of brainpathway_multisubject object\n');

    % Need to create a new object first with the same atlas
    bs2 = brainpathway_multisubject(bs.region_atlas);          % Construct a brainpathway object from an atlas object

    myfields = fieldnames(bs);  %{'voxel_dat' 'node_dat' 'region_dat' 'network_dat' 'partition_dat'
    myfields(strcmp(myfields, 'region_atlas')) = [];

    for i = 1:length(myfields)

        bs2.(myfields{i}) = bs.(myfields{i});
    end
    
    % NOTE: i do not understand why the below two lines are necessary! it
    % seems the contents of r and p don't get copied over, which is
    % puzzling. Below, commented out, is a clunky partially implemented
    % fix. Better to figure out the cause of the problem and then address,
    % rather than implementing cludgy solutions
    % -- Yoni Ashar, March 2020
    bs2.connectivity.regions.r = bs.connectivity.regions.r;
    bs2.connectivity.regions.p = bs.connectivity.regions.p;
%     
%         % recursive copy
%         if isstruct(bs.(myfields{i}))
%             mysubfields = fieldnames(bs.(myfields{i}));
%             for j = 1:length(mysubfields)
%                 b2.(myfields{i}).(mysubfields{j}) = bs.(myfields{i}).(mysubfields{j});
%             end
%             
%             if ~isempty(bs.(myfields{i})) & isstruct(bs.(myfields{i}).(mysubfields{j}))
%                 mysubsubfields = fieldnames(bs.(myfields{i}).(mysubfields{j}));
%                 for k = 1:length(mysubsubfields)
%                     b2.(myfields{i}).(mysubfields{j}).(mysubsubfields{k}) = bs.(myfields{i}).(mysubfields{j}).(mysubsubfields{k});
%                 end
%             end
%         end
%     end

    % select subjects
    bs2.subject_metadata(~wh_keep,:) = [];
    bs2.connectivity.regions.r(:,:,~wh_keep) = [];
    bs2.connectivity.regions.p(:,:,~wh_keep) = [];
    bs2.region_dat(~wh_keep) = [];
    
    for f = fieldnames(bs2.data_quality)'
        if ~isempty(bs2.data_quality.(f{1}))
            bs2.data_quality.(f{1})(~wh_keep) = [];
        end
    end
    if ~isempty(bs2.connectivity.nodes) && ~isempty(bs2.connectivity.nodes.r)
        bs2.connectivity.nodes.r(:,:,~wh_keep) = [];
        bs2.connectivity.nodes.p(:,:,~wh_keep) = [];
    end
    
    % TODO: expand subject selection to other fields, once they are getting
    % copied over correctly
end