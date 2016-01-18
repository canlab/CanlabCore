function [mysets,n_in_set,sets_by_vars,classes] = parcel_complete_sets(s,varargin)
% :Usage:
% ::
%
%     [mysets,n_in_set,sets_by_vars,classes] = parcel_complete_sets(s,['dounique','nofuzzy'])
%
% Find sets of elements in logical n x n matrix s in which all pairs in set are 'true' (1 in
% matrix s)
%
% Optional: input a distance/correlation/etc. matrix and threshold
%
% :Optional Inputs:
%
%   **'dounique':**
%        provides single-variable sets as well
%
%   **'nofuzzy':**
%        chooses closest set for each var so that each var
% can belong to only one set
%
%   **'threshold':**
%        followed by threshold thr for matrix s.
%
%        if thr is a string, thr will be evaluated on s: e.g., 's < 10'
%
%        in this case, thr should be a logical expr. involving s
%
%        if thr is a number, s >= thr will be evaluated
%
%        the resulting logical matrix will be used to determine sets
%
%   **'min' or 'max':**
%        works only with 'nofuzzy' option. 
%
%        if max, uses max to find most similar set; good if s is a
%        covariance matrix
%
%        if min, uses min to find closest set; good if s is a
%        distance matrix
%
%        default is max
%
% runs on Matlab 7.2 or higher %%%
%
% :Examples:
% ::
%
%    % Find sets of coordinates within 10 mm of one another
%    xyz = cat(1,cl{1}.mm_center);
%    d = pdist(xyz); d = squareform(d);
%    [mysets,n_in_set,sets_by_vars,classes] = parcel_complete_sets(d,'dounique','nofuzzy','threshold','s<10','min');
%
% ..
%    tor wager, nov 10, 2006
% ..

    keywords = varargin; % process inputs and initialize variables
    %keywords = {'dounique' 'nofuzzy' 'threshold'};
    classes = [];
    minmax = 'max';
    
    if any(strcmp(keywords,'min')), minmax = 'min'; end
    if any(strcmp(keywords,'max')), minmax = 'max'; end
    
    if any(strcmp(keywords,'threshold'))
        thr = keywords{find(strcmp(keywords,'threshold')) + 1}; % get threshold
        sorig = s;      % save original for 'nofuzzy' option (only needed here)

        if ischar(thr)
            eval(['s = ' thr ';']);
        else
            s = s >= thr;
        end
    else
        s = logical(s);
    end

    %  make sure diagonals are zeros, so that a var. can't correlate with itself
    s = s .* (1 - eye(size(s)));

    % this is a list of all the unique pairs of 'true' values
    [r,c] = find(triu(s));

    % number of pairs of 'true' values
    n = length(r);

    mysets = cell(1,n);

    n_in_set = zeros(1,n);
    sets_by_vars = false(n,size(s,2));

    % for each pair, get the set of variables connected with all members of the
    % set (recursive search)
    % get sets where all in set are 'related' (1 in matrix)
    % ---------------------------------------------------------------
    for i = 1:n
        mysets{i} = get_set;

        %     n_in_set(i) = length(mysets{i});

        sets_by_vars(i,mysets{i}) = 1;
    end

    % % [sets_by_vars,i] = unique(sets_by_vars,'rows');
    % % mysets = mysets(i);
    % % n_in_set = n_in_set(i);
    get_outputs_from_sets_by_vars;

    in_sets = sum(sets_by_vars,1);    % how many sets each var is in


    if any(strcmp(keywords,'dounique'))
        % for each original var, make sure it's in a set
        % if none, create unique set
        % ---------------------------------------------------------------
        newsets = find(in_sets == 0);
        nnew = length(newsets);
        [add_to_sets{1:nnew}] = deal(1);
        mysets = [mysets add_to_sets];
        n_in_set = [n_in_set ones(1,nnew)];
        nvars = size(s,1);
        add_to_sbyv = zeros(nnew,nvars);
        for i = 1:nnew
            add_to_sbyv(i,newsets(i)) = 1;
        end
        sets_by_vars = [sets_by_vars; add_to_sbyv];

    end

    if any(strcmp(keywords,'nofuzzy'))
        % for each original var, make sure it's in only one set
        % if > 1, choose closest set
        % ---------------------------------------------------------------
        manysets = find(in_sets > 1);

        for i = manysets
            select_closest_set;
        end

        % recompute outputs
        get_outputs_from_sets_by_vars;

        [classes,tmp] = find(sets_by_vars);

    end

    %%% End main function %%%
    %%% Matlab 7.2 or higher needed for nested functions, below %%%




    % ---------------------------------------------------------------
    % nested functions
    % ---------------------------------------------------------------

    function myset = get_set
        % start with the pair of correlations...
        % myset is the indices of variables that belong to a set.
        myset = [r(i) c(i)];
        r_this_set = s(myset,:);

        % ...and iteratively find others that correlate with all members of the
        % set.

        while true
            % add2set gives indices of variables that are correlated with all of the
            % others in the set
            add2set = find(all(r_this_set));
            %if any(add2set == 16), keyboard; end

            if isempty(add2set), break, end

            % add them one at a time to avoid adding incompatible vars
            myset = [myset add2set(1)];
            r_this_set = s(myset,:);
        end

    end

    function get_outputs_from_sets_by_vars

        [sets_by_vars] = unique(sets_by_vars,'rows');
        emptysets = sum(sets_by_vars,2) == 0;
        sets_by_vars(emptysets,:) = [];
        nsets = size(sets_by_vars,1);
        mysets = cell(1,nsets);

        for i = 1:nsets
            mysets{i} = find(sets_by_vars(i,:));
        end
        n_in_set = sum(sets_by_vars,2)';

    end

    function select_closest_set
        wh = find(sets_by_vars(:,i));   % the sets this var is in
        len = length(wh);
        meanset = zeros(1,len);

        % pick closest (farthest) candidate set in terms of sim.
        for j = 1:len                   % for each candidate set
            thisset = mysets{wh(j)};
            thisset(thisset == i) = []; % get members of candidate set not self

            % mean relationship values w/other members

            if exist('sorig','var')       % use actual values if we have them
                meanset(j) = mean(sorig(thisset,i));
            else
                meanset(j) = mean(s(thisset,i));
            end
        end

        % pick closest (farthest) candidate; in case of ties, pick first one
        switch minmax
            case 'max'
                [mymax,whmax] = max(meanset);
                myset = wh(whmax);
            case 'min'
                [mymax,whmax] = min(meanset);
                myset = wh(whmax);
        end

        % select vals in set; remove other memberships
        % remove from sets_by_vars
        sets_by_vars(:,i) = 0;
        sets_by_vars(myset,i) = 1;

    end

end  % end main function
