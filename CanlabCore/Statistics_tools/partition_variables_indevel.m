function best_classes = partition_variables_indevel(x, k, obsk)
% partition n x v matrix into k classes
%
% maximize within-class condition number or minimize cov. determinant
%
% minimize between-class condition number or maximize cov. determinant...
%
% Simulate two-class data
% ::
%
%    nvars = 100; nsubj = 20; corval = .5;
%    S = eye(nvars./2) + corval*(1 - eye(nvars./2));
%    S = blkdiag(S,S); S(S==0) = corval;
%    x = mvnrnd(zeros(1,nvars), S, nsubj); det((corrcoef(x)))
%    figure; imagesc(corrcoef(x)); colorbar



    if nargin < 3, obsk = 1; end  % classes of observations
    v = size(x,2);
    npossible = k ^ v;

% create population of possible solutions
% k^v possible solutions, where v is variables and k is classes
% we can only do this with very small sets of variables

    
    if npossible <= 10 ^ 5
        fprintf(1,'Could do exhaustive search\n');
        clist = repmat({1:k}, 1, v);
        pop = combvec(clist{:});
    end

    % ------------------------------------------------------------
    % set up GA
    % ------------------------------------------------------------
    gensize = 1000;
    numgen = 30;
    fprintf(1,'Running GA with %3.0f organisms and %3.0f generations.\n', gensize, numgen);

    % this defines the fitness function to maximize
    % the variable input parameter is class assignment
    % other parameters (data and max # of classes) are fixed
    if obsk == 1
        objective_function = @(class) mean_correl_fitness(class, x, k);
    else
        objective_function = @(class) twoway_correl_fitness(class, x, k);
    end
    
    inputs{1} = ones(v, 1);
    inputs{1}(1:k) = (1:k)';

    % ------------------------------------------------------------
    % run GA
    % ------------------------------------------------------------

    [best_classes,fit,beff,in] = tor_ga(gensize,numgen,inputs,objective_function,'discrete');
    best_classes = best_classes{1};

    if obsk > 1

        n = size(x,1);  % # observations

        best_obsclasses = best_classes(1:n);
        best_classes = best_classes(n+1:end);

    end

    % ------------------------------------------------------------
    % plots, etc.
    % ------------------------------------------------------------
    %%%fit = determinant_fitness(class, x, k)
    [bestsort, wh] = sort(best_classes);
    xsort = x(:, wh);

    classavg = [];
    for i = unique(best_classes)'
        classavg(:,end+1) = mean(x(:, best_classes == i),2);
    end

    figure; subplot(3,2,1)
    imagesc(xsort); xlabel('Variables'); ylabel('Cases');
    subplot(3,2,2);
    imagesc(classavg); xlabel('Class'); ylabel('Cases');
    subplot(3,2,3);
    imagesc(bestsort'); xlabel('Class');

    subplot(3,2,5);
    imagesc(corrcoef(xsort)); xlabel('Variables'); ylabel('Variables');



    figure; subplot(3,2,1)
    imagesc(x); xlabel('Variables'); ylabel('Cases');
    subplot(3,2,2);
    imagesc(classavg); xlabel('Class'); ylabel('Cases');
    subplot(3,2,3);
    imagesc(best_classes'); xlabel('Class');

    subplot(3,2,5);
    imagesc(corrcoef(x)); xlabel('Variables'); ylabel('Variables');


end






function fit = determinant_fitness(class, x, k)
    % higher is better


    c = corrcoef(x);

    bindx = 1;

    % get non-empty classes; these are omitted from wiclass and btclass
    % they do not influence fitness score at all

    valid_classes = unique(class)';
    n_valid_classes = length(valid_classes);

    if iscol(valid_classes), valid_classes = valid_classes'; end

    % initialize this to ones so empty classes contribute -log(1) = 0 to
    % overall fitness
    wiclass = ones(1, k);

    btclass = [];

    for i = valid_classes  % for each class
        wh = class == i;

        wiclass(i) = det(c(wh,wh));   % determinant of within-class vars; minimize this.  range: of det: btwn 0 and 1

        for j = i+1 : k

            if any(valid_classes == j)
                % if this class is valid
                whj = class == j;

                btwncorr = abs(c(wh,whj));

                btclass(bindx) = -mean(log(btwncorr(:)));  % product of between-class vars; minimize this; maximize -log
                % don't know how to scale this relative to wiclass; range of mean: btwn 0 and 1

                bindx = bindx + 1;

            else
                % skip
            end
        end

    end

    wiclass = -log(wiclass);

    nbtwn_els = (bindx - 1);
    fit = sum([wiclass ./ n_valid_classes   btclass ./ nbtwn_els]);

end





function fit = twoway_correl_fitness(class, x, k)
    % class is concatenated: observation class, then variable class
    
    v = size(x,2); % # variables
   
    n = length(class) - v;  % # observations
    
    obsclass = class(1:n);
    class = class(n+1:end);
    
    valid_classes = unique(obsclass)';
    
    for i = valid_classes
        
        whobs = obsclass == i;  % which observations for this obs. grouping
        fit(i) = mean_correl_fitness(class, x(whobs, :), k);
        
        npooled(i) = sum(whobs);
    end
    
    fit = sum(fit .* npooled) ./ sum(npooled);
    
end
    
    




function fit = mean_correl_fitness(class, x, k)
    % higher is better

    c = corrcoef(x);

    % ------------------------------------------------------------
    % setup
    % ------------------------------------------------------------
    
    bindx = 1;

    % get non-empty classes; these are omitted from wiclass and btclass
    % they do not influence fitness score at all

    valid_classes = unique(class)';
    n_valid_classes = length(valid_classes);

    if iscol(valid_classes), valid_classes = valid_classes'; end

    % initialize this to zeros so empty classes 0 to
    % overall fitness
    wiclass = zeros(1, k);
    nwi = zeros(1, k);

    btclass = [];

    % ------------------------------------------------------------
    % loop through classes
    % ------------------------------------------------------------
    
    for i = valid_classes  % for each class
    
    % within-class
    % ------------------------------------------------------------
    wh = class == i;

        wicorr = c(wh,wh);
        nwi(i) = sum(wh);

        if nwi(i) == 1
            wiclass(i) = 1;
        else
            wicorr = wicorr - eye(size(wicorr));
            wicorr = squareform(wicorr);

            wiclass(i) = mean(wicorr);   % correl of within-class vars; maximize this.
        end
        
        for j = i+1 : k
            % between-class
            % ------------------------------------------------------------

            if any(valid_classes == j)
                % if this class is valid
                whj = class == j;

                btwncorr = c(wh,whj);
                btwncorr = btwncorr(:);

                btclass(bindx) = mean(btwncorr);  % condition # of btwn-class vars; minimize this.
                nbt(bindx) = length(btwncorr);

                bindx = bindx + 1;

            else
                % skip
            end
        end

    end

    % finish up
    % ------------------------------------------------------------

    nbtwn_els = (bindx - 1);

    wifit = sum( wiclass .* nwi ) ./ sum(nwi);  % pos is good, zero or neg is bad , range -1 - 1
    %btfit = sum( abs(btclass) .* nbt ) ./ sum(nbt);   % zero is good, pos or neg is bad , range 0 - 1
    btfit = sum( (btclass) .* nbt ) ./ sum(nbt);   % neg better than 0 better than pos

    % total range: -2 to 1; random start should be around 0
    fit = wifit - btfit;

end
