function newcl = merge_nearby_clusters(cl, thr, varargin)
% Merge sets of clusters whose centers are all within thr mm of each other
% uses parcel_complete_sets.m
%
% :Usage:
% ::
%
%    newcl = merge_nearby_clusters(cl, thr)
%
% :Example:
% ::
%
%    % The command below runs the function recursively until all clusters are
%    % greater than thr mm apart
%    newcl = merge_nearby_clusters(cl, thr, 'recursive')
%
% ..
%    tor wager, nov 06
% ..
%


    if isempty(cl), newcl = []; return, end

    dorecursive = 0;
    % run this function recursively if specified
    if any(strcmp(varargin, 'recursive'))
        d = 0;
        while any(d < thr)
            cl = merge_nearby_clusters(cl, thr);
            d = pdist(cat(1, cl.mm_center));
        end
        
        newcl = cl;

        return
    end


    evalstr = ['s < ' num2str(thr)];


    % get distances based on cluster centers
    xyz = cat(1, cl.mm_center);
    d = pdist(xyz);
    d = squareform(d);

    [mysets, n_in_set, sets_by_vars, classes] = parcel_complete_sets(d, 'dounique', 'nofuzzy', 'threshold', evalstr, 'min');

    for i = 1:length(mysets)
        wh = mysets{i}; % which clusters to merge into new cluster i

        % start the new cluster
        newcl(i) = cl(wh(1));

        % get target sizes to match
        % updates N (names), nfields
        get_fields_to_cat;

        for j = 2:length(wh)
            for k = 1:nfields
                [dat, sz] = get_sz_and_data(cl(wh(j)) , N{k});

                update_newcl;

            end

            update_titles;
        end

        newcl(i).numVox = size(newcl(i).XYZ, 2);
        newcl(i).mm_center = center_of_mass(newcl(i).XYZmm, newcl(i).Z);
    end

    % NESTED
    %
    %

    function get_fields_to_cat
        N = fieldnames(cl);
        nfields = length(N);
        for k = 1:nfields
            [dat, sz(k, :)] = get_sz_and_data(newcl(i), N{k});
        end
        wh_horzcat = find(sz(:, 2) == newcl(i).numVox);
        N = N(wh_horzcat);

        % remove specific problem fields
        whr = find(strcmp(N, 'M')); if ~isempty(whr), N(whr) = []; end
        whr = find(strcmp(N, 'mat')); if ~isempty(whr), N(whr) = []; end
        whr = find(strcmp(N, 'dim')); if ~isempty(whr), N(whr) = []; end
        whr = find(strcmp(N, 'voxSize')); if ~isempty(whr), N(whr) = []; end

        nfields = length(N);
    end

    function update_newcl
        if sz(1) == size(newcl(i).(N{k}), 1)
            newcl(i).(N{k}) = [newcl(i).(N{k}) dat];
        else
            warning(['Size mismatch in field ' N{k}])
        end
    end

    function update_titles
        if isfield(newcl, 'title'), newcl(i).title = str2mat(newcl(i).title, cl(wh(j)).title); end
        if isfield(newcl, 'name'), newcl(i).name = str2mat(newcl(i).name, cl(wh(j)).name); end
        if isfield(newcl, 'shorttitle'), newcl(i).shorttitle = str2mat(newcl(i).shorttitle, cl(wh(j)).shorttitle); end
    end

end


function [dat, sz] = get_sz_and_data(cl, fldname)
    dat = cl.(fldname); % data for this field
    sz = size(dat);         % size of field data
end
