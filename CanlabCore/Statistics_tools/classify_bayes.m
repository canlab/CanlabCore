function [corrclass, taskclass, realclass, likeratio, m, misclass,ptask,indx,cl] = classify_bayes(meth,y,Xi,varargin)
% :Usage:
% ::
%
%     [corrclass, taskclass, realclass, likeratio, m, misclass,ptask,indx,cl] = classify_bayes(meth,y,Xi,[var args])
%
% :Inputs:
%
%   **y:**
%        observations (e.g., studies) x variables (e.g., brain voxels)
%
%        y = SOMResults.dat';
%
%   **Xi:**
%        Xi = MC_Setup.Xi(:,1:2);
%
%   **meth:**
%        can be: {'linear','diagLinear','quadratic','diagQuadratic','mahalanobis'}
%
%        discriminant analysis
%        'bayes' : simple Bayes posterior prob classifier
%
% :Optional inputs:
%
% Feature selection parameters:
%
%   **selectivity_cutoff:**
%        max probability of task given a response in a variable,
%        divided by number of tasks.
%
%        1 = variable must exceed .5 for 2 tasks, .2 for 5 tasks, etc.
%
%        1.5 = .3 for 5 tasks, .75 for 2 tasks, etc.
%
%        0 = no selectivity
%
%   **activation_cutoff:**
%        max proportion of studies of some type that produced a
%        response in a variable (e.g., voxel)
%
%        .1 is default
%
%        0 is no selectivity
%
% :Examples:
% ::
%
%    [corrclass, taskclass, realclass, likeratio, m, misclass] = ...
%    classify_bayes('bayes',MC_Setup.unweighted_study_data',Xi,'selectivity_cutoff',1); corrclass
%
%    [corrclass, taskclass, realclass, likeratio, m, misclass] =
%    classify_bayes('bayes',MC_Setup.unweighted_study_data',Xi,'selectivity_cutoff',1,'activation_cutoff',.08); corrclass
%
% Batch modes:
% Permutation test: input 'permtest' followed by number of permutations
% Do not request more than one output variable
% [corrclass_nullhyp] = classify_bayes('bayes',dat,Xi,'selectivity_cutoff',1.8,'activation_cutoff',.02,'permtest',5);
%
% To get mask index of which areas meet feature selection:
% whsave = sum(indx > 0, 2) > 0;
%
% ..
%    tor wager, 10/2/06
% ..

    % select features first, instead of in x-validation (makes x-val
    % invalid; don't do it)
    selectfirst = 0;

    dofeature = 1;
    selectivity_cutoff = 1.5;
    activation_cutoff = .1;
    dopca = 0;
    doxval = 1;
    doplot = 1;

    for i = 1:length(varargin)
        if isstr(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'nofeature', dofeature = 0;
                case 'reduction', dopca = 0;
                case 'noxval', doxval = 0;
                case 'noplot', doplot = 0;
                case 'selectfirst', selectfirst = 1;
    
                    % functional commands
                case 'selectivity_cutoff', selectivity_cutoff = varargin{i+1};
                case 'activation_cutoff', activation_cutoff = varargin{i+1};

                    % batch modes
                case 'permtest'  % batch permutation test on shuffled task indicators
                    nperms = varargin{i+1};
                    if ~exist('varargin','var'), varargin = {}; end
                    corrclass = batch_permtest(nperms,meth,y,Xi,varargin);

                    return

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    origvars = size(y,2);
    [nstudies,ntasks] = size(Xi);
    realclass = indic2condf(Xi);

    % eliminate no-variance voxels
    % -----------------------------------------------------------------
    [y,whsave] = eliminate_constant(y);

    nvars = size(y,2);
    nalltasks = sum(sum(Xi));

    fprintf(1,'\n')
    str = ('Inital computation'); disp(str)

    % get max and min ptask, to avoid taking log of 1 or 0
    % -----------------------------------------------------------------
    % can't take log of 0, so shrink these extreme values
    sy = sum(y);
    ptask = (Xi' * y) ./ repmat(sy,ntasks,1);
    mymax = max(ptask(ptask < 1));
    mymin = min(ptask(ptask > 0));

    erase_string(str)

    % y is all data
    % trainset = training set, studies x voxels (regions)
    % testvec = test vector, 1 x voxels (regions)

    selectivity_cutoff = selectivity_cutoff./ntasks;
    str = sprintf('Feature selection: selectivity (/ntasks): %3.2f, activation: %3.2f\n',selectivity_cutoff,activation_cutoff);
    disp(str)

    if dofeature && selectfirst
        [y,whsave] = select_features(ptask,Xi,y,selectivity_cutoff,activation_cutoff,whsave);
        nvars = size(y,2);
    end

    if dopca
        % reduce data
        % -----------------------------------------------------------------
        y = pca_reduction(y);
    end

    maxlike = zeros(nstudies,1);
    taskclass = zeros(nstudies,1);
    likeratio = zeros(nstudies,1);

    if doxval
        % -----------------------------------------------------------------
        fprintf(1,'Cross-validating: 00000');
        for i = 1:nstudies

            if mod(i,10) == 0, fprintf(1,'\b\b\b\b\b%05d',i); end
            testvec = y(i,:)';
            trainset = y;
            trainset(i,:) = [];

            X = Xi;
            X(i,:) = [];

            sumtasks = sum(trainset); sumtasks(sumtasks==0) = 1;  % avoid warning for empty tasks
            ptask = (X' * trainset) ./ repmat(sumtasks,ntasks,1);
            % feature selection
            if dofeature && ~selectfirst
                trainset = select_features(ptask,X,trainset,selectivity_cutoff,activation_cutoff,whsave);
            end

            [taskclass(i),maxlike(i),likeratio(i)] = get_class(meth,X,trainset,testvec,ntasks,mymax,mymin,ptask);
        end

    else
        % -----------------------------------------------------------------
        % no x validation
        ptask(ptask == 1) = mymax;
        ptask(ptask == 0) = mymin;

        for i = 1:nstudies
            testvec = y(i,:)';
            [taskclass(i),maxlike(i),likeratio(i)] = choose_most_likely(ptask,testvec); % taking the log inside here is not efficient
        end
    end %doxval



    % finish up: final features
    sumtasks = sum(y); sumtasks(sumtasks==0) = 1;  % avoid warning for empty tasks
    ptask = (Xi' * y) ./ repmat(sumtasks,ntasks,1);
    if dofeature
        [y,whsave,nvars] = select_features(ptask,Xi,y,selectivity_cutoff,activation_cutoff,whsave);
    end

    sumtasks = sum(y); sumtasks(sumtasks==0) = 1;  % avoid warning for empty tasks
    ptask = (Xi' * y) ./ repmat(sumtasks,ntasks,1);

    [m,dp,corr,far,misclass] = confusion_matrix(realclass,taskclass);

    wh = (realclass ~= 0);
    n = sum(wh);
    corrclass = 1 - sum(misclass) ./ n;

    % get class vectors for output, in original image length
    [mx,r] = max(ptask,[],1);

    indx = zeros(origvars,ntasks);
    for i = 1:ntasks
        thistask = (r == i); % & (max(ptask) > ptaskcutoff);
        indx(whsave(thistask),i) = ptask(i,thistask)';
    end

    cl = [];
    if doplot
        cl = classify_viz_regions(indx);
    end

    return


    % notes


    %%%no  prior_task = sum(Xi) ./ sum(Xi(:));
    %ptask = ptask .* repmat(prior_task',1,nvars);

    % this is exactly proportional to above but wrong.
    %pa = repmat(sum(y) ./ nstudies,ntasks,1);
    %ptask = (Xi' * y) ./ nalltasks ./ pa;



    % ====================================================================
    % --------------------------------------------------------------------
    %
    % Support functions
    %
    % --------------------------------------------------------------------
    % ====================================================================

function [taskclass,maxlike,likeratio] = get_class(meth,X,trainset,testvec,ntasks,mymax,mymin,varargin)
    switch meth
        case 'bayes'
            % proportion of activations that came from each task
            % p task | activation, with flat priors for task
            % posterior probability of task | activation in voxel (Region)
            if length(varargin) > 0
                % do this for x-validation; we already computed ptask
                ptask = varargin{1};
            else
                sumtasks = sum(trainset); sumtasks(sumtasks==0) = 1;  % avoid warning for empty tasks
                ptask = (X' * trainset) ./ repmat(sumtasks,ntasks,1);
            end

            %(xy / sumx) * (sumx./allsum) ./ (sumy./nstudies)
            % add task priors
            %prior_task = sum(X) ./ sum(X(:));
            %ptask = ptask .* repmat(prior_task',1,nvars);

            % gives posterior probability
            % %                 pa = repmat(sum(trainset) ./ (nstudies-1),ntasks,1);
            % %                 ptask = (Xi' * trainset) ./ (nalltasks-1) ./ pa(j);

            ptask(ptask == 1) = mymax;
            ptask(ptask == 0) = mymin;

            % choose most likely task
            [taskclass,maxlike,likeratio] = choose_most_likely(ptask,testvec);

% %         case 'naivebayes'
% %             [taskclass log_joint best_log map] = classify_naive_bayes('test', Y, bayes_model);
            

        case {'linear','diagLinear','quadratic','diagQuadratic','mahalanobis'}
            condf = realclass;
            condf(i) = [];
            %[taskclass(i),err,ptask,logp] = classify(testvec',trainset,condf,meth);
            [taskclass(i),err] = classify(testvec',trainset,condf,meth);

        otherwise error('Unknown method.')
    end

    return


function [taskclass,maxlike,likeratio] = choose_most_likely(ptask,testvec)


    ptask = log(ptask);    % likelihoods; summing these = taking product of probabilities

    taskprob = ptask * testvec;

    % choose most likely task
    [maxlike,taskclass] = max(taskprob);
    if all(taskprob == 0)
        likeratio = NaN;
    else
        likeratio = max(taskprob) ./ min(taskprob);
    end

    return


function [y,wh] = eliminate_constant(y)
    str = sprintf('Eliminating constant variables'); disp(str);
    %whomit = find(all(y == repmat(mean(y),nstudies,1)));
    tmp = sum(y);
    whomsave = ~(tmp == 0); %| tmp == nstudies);

    y = y(:,whomsave);
    wh = find(whomsave);
    erase_string(str);

    return



function [y,whsave,nvars] = select_features(ptask,Xi,y,selectivity_cutoff,activation_cutoff,whsave)
    % select features
    % -----------------------------------------------------------------

    maxp = max(ptask);  % max prob of task across vars
    nvars = size(y,2);
    pt = (Xi' * y) ./ repmat(sum(Xi)',1,nvars);  % proportion of each task that activated in an area
    ptmax = max(pt);
    whomsave = maxp > (selectivity_cutoff) & ptmax > activation_cutoff;

    if sum(whomsave) == 0
        disp('Warning: No variables meet feature selection criteria.');
        whomsave = (maxp == max(maxp));
    end

    whsave = whsave(whomsave);  % overall indices of saved voxels
    y = y(:,whomsave);
    ptask = ptask(:,whomsave);
    nvars = size(y,2);
    %fprintf(1,'Selected: %3.0f ',nvars);
    % get weights on voxels based on n activations.... ******

    return


function score = pca_reduction(x)

    str = sprintf('\nData reduction'); disp(str);
    [pc,score,latent] = princomp(full(x),'econ');

    % %         if issparse(x)
    % %             [score,S,V] = svds(x,'econ');
    % %         else
    % %             [score,S,V] = svd(x,'econ');
    % %         end
    % %         latent = diag(S);

    % Eigenvalue plot
    figure('Color','w'), set(gca,'FontSize',18),bar(latent),
    xlabel('Components'),ylabel('Variance')

    erase_string(str);
    num_to_save = sum(latent>1);
    num_to_save = input(sprintf('Enter num. to save (%3.0f are > 1) : ',num_to_save));

    score = score(:,1:num_to_save);

    return



    % ====================================================================
    % --------------------------------------------------------------------
    %
    % Batch mode functions: permutation test
    %
    % --------------------------------------------------------------------
    % ====================================================================

function corrclass = batch_permtest(nperms,meth,y,Xi,varinputs)

    % build f call based on inputs
    str = 'ccperm = classify_bayes(meth,y,Xiperm,''noplot''';
    if ~isempty(varinputs)
        for i = 1:length(varinputs)
            % remove call to permtest
            toadd = varinputs{i};
            if ischar(toadd) && strcmp(toadd,'permtest')
                varinputs{i} = []; varinputs{i+1} = [];
            end

            % add input
            toadd = varinputs{i};

            if isempty(toadd)
                % nothing to do.
            elseif ischar(toadd)
                toadd = ['''' toadd ''''];
            else
                toadd = num2str(toadd);
            end

            if isempty(toadd)
                % nothing to do.
            else
                str = [str ',' toadd];
            end
        end
    end
    str = [str ');'];

    fprintf(1,'Permutation test on permuted Xi.\nPermutations: %3.0f\n',nperms);
    fprintf(1,'Function call:\n%s\n',str);

    corrclass = zeros(1,nperms);
    for i = 1:nperms
        Xiperm = shuffles(Xi);
        eval(str);
        corrclass(i) = ccperm;
        
        if mod(i,10) == 0
            disp('Saving results so far in tmp_corrclass_permuted');
            save tmp_corrclass_permuted corrclass
        end
    end

    return



    % ====================================================================
    % --------------------------------------------------------------------
    %
    % Other utility functions
    %
    % --------------------------------------------------------------------
    % ====================================================================

function erase_string(str)
    fprintf(1,repmat('\b',1,length(str)+1)); % erase string
    %fprintf(1,'\n');
    return
