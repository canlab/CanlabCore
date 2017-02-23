function rankout = rankdata(m, varargin)
    % ranks = rankdata(m)                % default: midrank method
    % ranks = rankdata(m, 'nomidrank')   % do not use midrank
    %
    % ranks a vector, returns ranks in place of values
    %
    % Tor Wager
    %
    % Modified April 2007 to use the midrank method to handle ties.
    % Modified May 2007 to allow option to rank unique values (useful in some
    % functions, such as computing Kendall's tau); see correlation.m



    [nobs, ncols] = size(m);

    rankout = zeros(nobs, ncols);

    for i = 1:ncols

        rankout(:, i) = getranks(m(:, i), nobs, varargin{:});

    end

end   % end main function






function ranks = getranks(m, nobs, varargin)

    ranks = zeros(nobs, 1);

    % first col of a is data, 2nd is indices
    % c is sorted ranks
    % ranks is ranks in original format

    [a,b]=sortrows([m (1:length(m))'],1);

    c = (1:length(m))'; ranks(a(:,2)) = c;



    % handle ties

    if ~isempty(varargin) && strcmp(varargin, 'nomidrank')
        % ranks on unique values

        m2 = unique(m);
        [a,b]=sortrows([m2 (1:length(m2))'],1);

        c = (1:length(m2))'; ranks1(a(:,2)) = c;

        for i = 1:length(m2)
            ranks(m == m2(i)) = ranks1(i);
        end


    else
        % midrank method

        uniqm = unique(m);

        if length(uniqm) ~= nobs

            for i = 1:length(uniqm)

                wh = a(:,1) == uniqm(i);

                if sum(wh) > 1

                    c(wh, 1) = mean(c(wh, 1));

                end
            end

            ranks(a(:,2)) = c;

        end

    end

end

