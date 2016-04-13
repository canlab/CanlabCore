function [best_params,fit,beff,in,isconverged] = tor_ga(gensize,numgen,inputs,ofun,varargin)
% :Usage:
% ::
%
%     [best_params,fit,beff,in,isconverged] = tor_ga(gensize,numgen,inputs,ofun,[optional in any order: genfun,fixed inputs,cmd strings])
%
% :Inputs:
%
%   a cell array describing the inputs to the optimization function
%   (parameters to be optimized).
%
% Each cell of inputs is a p x q matrix of parameters.
% p and q are arbitrary, as each organism is described by a p x q
% matrix...but the objective function must be able to handle inputs in
% the format you provide.
% Internally, a set of 'organisms' is created that is
% params x params x organisms (3-D).
% This matrix is subject to crossover across orgs. separately for each
% cell.
% If inputs is a p x q matrix, it will be placed in a single cell.
%
% By default, it is only necessary to enter a single set of params
% for an example organism.  The range of those input values is used to
% generate random starting values for each organism.
%
% There can be more than one set of
% parameters that are combined in some way by ofun to produce a fitness
% value.  if there is more than one set of input parameters,
% inputs should be entered as a cell array, one cell per input.
% inputs should be in ORDER of inputs entered to ofun!
%
% **RECOMMENDED**
%
% If you enter each cell of inputs as a 3-D array so that inputs(:,:,1)
% is the min acceptable value for each param and inputs(:,:,2) is the
% max acceptable value, then the ga will create a series of organisms
% at start that evenly span the range of the multivariate parameter
% space, with the spacing between values determined by the gensize.
% This can provide a huge advantage in efficiency for the GA.
% This is the idea behind the "Sobol sequence," which chooses values
% that evenly span a multivariate space.
% With this option, if gensize is sufficiently large and the param
% space is sufficiently small, then the ga may find the correct
% solution on the first iteration.
% However, it is not likely to work well if the num. params >> gensize
%
% :Examples:
% ::
%
%    start = [-15 -15];
%    start(:,:,2) = [15 15];
%    [best_params,fit,beff,in] = tor_ga(324,30,{start},objfun_ga,'genconverge',5);
%
% **ofun**
%   - the objective function that combines the inputs.
%
% There are two options for passing this in:
%   1) enter the name of the function as a string.  the program creates a handle for the
%      function, and evaluates it using inputs specified in the inputs variable.
%      In this case, pass in fixed inputs after ofun, in the varargin fields
%      fixed inputs  optional, fixed inputs that do not change!  same structure
%      as inputs.
%
%   2) You can also enter ofun as a function handle durectly, with fixed inputs already
%      embedded before running the program.
%      The function should take as input a param list, and return fitness.
%      e.g.,
%      ::
%
%         objhan = @(params) my_function_name(params,fixed_inputs1,fixed_inputs1,fixed_inputs1);
%         objhan = @(wh) prospect_organism(ceil(wh),pop,truep,iter);
%
% Pass in objhan as the 'ofun' input argument
%
% **genfun**
%   - [optional] input param generation function
%
% A function handle that generates a parameter set for each organism
%
% **command strings**
%   - 'noverbose'   turn off verbose reporting and plots
%   - 'genconverge' followed by integer x: converge if no change in last x generations
%
% ..
%    by Tor Wager, Last updated: Feb 2007, then June 2010 to add seeds
% ..
%
% :Examples:
%
% Example for fitting indscal model:
% ::
%
%    inputs{1} = X; fixin{1} = sp; fixin{2} = B1; fixin{3} = B2;
%    tor_ga(30,10,inputs,'indscalf',fixin);
%
%    W = rand(size(W)); W(1,:) = [10 10];, W(2,:) = [-10 -10];
%    inputs{2} = W;
%
% Example: Optimize gambles for prospect theory model
% See prospect_optimize_design.m for definition of population of
% gambles from which to draw (pop), truep, iter (all fixed inputs)
% ::
%
%    objhan = @(wh) prospect_organism(ceil(wh),pop,truep,iter);
%    genfun = @() randsample(gindx,ntrials,'true')';
%    [best_params,fit,beff,in] = tor_ga(5,3,wh,objhan,genfun);
%
% Using string inputs to control behavior:
% ::
%
%    [best_params,fit,beff,in] = tor_ga(300,30,{[15; -15]},objfun_ga,'genconverge',5,'noverbose');

    t0 = clock;

    % --------------------------------------------------------------------
    % * set up inputs
    % --------------------------------------------------------------------

    % objective function
    % two modes: 1)  given string, construct string to evaluate with fixed inputs
    %            2)  given function handle with fixed inputs embedded, evaluate
    %            directly
    if ischar(ofun)
        evalmode = 'string';
        eval(['fun = @' ofun ';']), disp(['Objective function is ' ofun ])
    else
        evalmode = 'handle';
        fun = ofun;

    end

    % format inputs correctly
    if ~iscell(inputs), tmp=inputs; clear inputs; inputs{1} = tmp; end

    switch evalmode
        case 'string'
            % inputs - create string that tells feval what arguments to put in


            estr = 'f = feval(fun,in{1}(:,:,j)';
            for i = 2:length(inputs)
                estr = [estr ',in{' num2str(i) '}(:,:,j)'];
            end
        case 'handle'
    end

    % optional inputs - fixed, non-optimized inputs
    % --------------------------------------------------------------------

    doverbose = 1;
    genconverge = 20;
    paramtype = 'continuous';
    dostochastic = 0;
    seeds = {};

    if ~isempty(varargin)
        for i = 1:length(varargin)

            if isempty(varargin{i})
                % do nothing

            elseif strcmp(class(varargin{i}),'function_handle')
                % this is the org. generation function
                genfun = varargin{i};

            elseif ischar(varargin{i})
                % command string
                switch varargin{i}
                    case 'noverbose', doverbose = 0;
                    case 'genconverge', genconverge = varargin{i+1}; varargin{i+1} = [];
                        
                    case {'integer', 'discrete'}, paramtype = 'discrete';
                        
                    case 'stochastic', dostochastic = 1;
                        
                    case {'seed', 'seeds'}, seeds = varargin{i+1};    varargin{i+1} = [];
                        
                    otherwise
                        disp('Warning: Unknown string input.')
                end

            else
                % this is fixed inputs
                fixin = varargin{i};
                if ~iscell(fixin), tmp=fixin; clear fixin; fixin{1} = tmp; end
                
                for j = 1:length(fixin)
                    estr = [estr ',fixin{' num2str(j) '}'];
                end

            end

        end
    end

% Check
if strcmp(paramtype, 'discrete') && dostochastic
    error('Parameter type must be continuous in order to use stochastic noise option.');
end

    
    if doverbose

        fprintf(1,'___________________________________________________________\n')
        switch evalmode
            case 'string'
                estr = [estr ');'];
                disp(['Evaluation string is ' estr])
            case 'handle'
                disp('Will evaluate this objective function: ')
                disp(fun)
        end
        
        disp(['GA thinks parameters are ' paramtype])
        ynstr = {'No' 'Yes'};
        fprintf('Status of stochastic noise addition on each generation: %s\n', ynstr{dostochastic + 1});
        fprintf('Will consider GA converged after %3.0f generations with no change\n', genconverge);
        
        fprintf(1,'___________________________________________________________\n')
    end

    % --------------------------------------------------------------------
    % * create start state
    % --------------------------------------------------------------------
    % create start state for fixed inputs

    for i = 1:length(inputs)

        if exist('genfun','var')
            % we have a custom organism-generation function handle
            if doverbose, disp('Generating starting values with custom (user input) organism generation function.'); end
            
            for j = 1:gensize
                in{i}(:,:,j) = genfun();
            end
        elseif size(inputs{i},3) > 1
            % we have a range of min/max values entered.  create param
            % values that span range.
            if doverbose, disp('Creating starting values that evenly span parameter space.'); end
            
            [in{i},gensize] = sobol_start_state(inputs{i},gensize,doverbose);
                        
            % we need these to limit noise later
            mininputs = repmat(inputs{1}(:,:,1), [1 1 gensize]);
            maxinputs = repmat(inputs{1}(:,:,2), [1 1 gensize]);

            
        else
            % we have only single starting value; use total range to create
            % random start state
            if doverbose, disp('Generating starting values with random parameter values.'); end

            in{i} = random_start_state(inputs{i},gensize,paramtype);
            
          

        end

    end
    
    % add custom seeds - tor added June 2010
    % --------------------------------------------------------------------
    if exist('seeds', 'var') && ~isempty(seeds)

        for i = 1:length(seeds)
            in{1}(:, :, end+1) = seeds{i};
        end

        gensize = size(in{1}, 3);

    end

    % --------------------------------------------------------------------
    % * set up info for adding randomness in case we get stuck later
    % --------------------------------------------------------------------
    % matrix of multiplier values for each param in first cell to add noise
    % (for continuous parameter spaces only)
    % or some needed info for discrete re-randomization
    if strcmp(paramtype, 'continuous')
        r = range(cat(3,in{:}),3); noisestd = .1 .* r;
        noisestd = repmat(noisestd,[1 1 size(in{1},3)]);
        
    elseif strcmp(paramtype, 'discrete')
        [m, n] = size(inputs{1});
        possiblevals = unique(inputs{1});  % should enter all possible values in input matrix
        nvals = length(possiblevals);
    end
    
    % --------------------------------------------------------------------
    % * iterate
    % --------------------------------------------------------------------

    if doverbose
        f1 = create_figure('Fitness by generation'); set(gcf,'Color','w'); hold on; title('Max fitness'),xlabel('Generation')
        
        f2 = create_figure('Params by generation'); set(gcf,'Color','w'); hold on; title('Parameter estimates'),xlabel('Generation')
        
        tmp = in{1}(:,:,1); tmp = tmp(:);
        all_bp = NaN .* zeros(numgen, length(tmp));
    end

    beff = NaN .* zeros(1,numgen);
    gentime = NaN .* zeros(1,numgen);
    isconverged = 0;

    fit = NaN .* zeros(gensize,numgen);
    

        

    for i = 1:numgen

        t1 = clock;

        if doverbose
            str = sprintf('Generation: %3.0f  ',i); fprintf(1,str);
            str = sprintf('Eval. fitness '); fprintf(1,str);
        end

        for j = 1:gensize

            % --------------------------------------------------------------------
            % * make models
            % --------------------------------------------------------------------

            % --------------------------------------------------------------------
            % * test models
            % --------------------------------------------------------------------
            switch evalmode
                case 'string'
                    % if using function eval string
                    eval(estr)          % e.g., f = fun(in{1}(:,:,j),fixin{1},fixin{2},fixin{3});
                    fit(j,i) = f;       % fitness of each model in each generation
                case 'handle'
                    fit(j,i) = fun(in{1}(:,:,j)); % in is matrix of params for all organisms, with 3rd-D, j, indexing org
            end


        end

        if doverbose
            erase_string(str);
            str = sprintf('Crossover '); fprintf(1,str);
        end

        % --------------------------------------------------------------------
        % * save best and crossover
        % --------------------------------------------------------------------
        eff = fit(:,i)';

        % run crossover for each separate input matrix to be recombined
        % best_params contains best parameters of this generation
        
        for j = 1:length(inputs)
            % replace in with crossover; b = index of best before xover
            % best_params is best of this generation
            % index of best after xover is 1 (hard-coded)
            % dostochastic indicates whether (some) children have some random
            % noise added to parent values
            
            % Note: inputs are shuffled after this point, and eff no longer
            % applies.  best input is saved in position 1, and in
            % best_params{j}
            [in{j}, b, best_params{j}] = xover(in{j},eff);
        end

        % --------------------------------------------------------------------
        % * add noise to perturb if needed 
        % --------------------------------------------------------------------
        
        if dostochastic
            % We've asked to add random noise to new organisms on every
            % generation.
            first_rand = ceil(gensize ./ 2);
            in{1}(:,:,first_rand:end) = in{1}(:,:,first_rand:end) + noisestd(:,:,first_rand:end) .* randn(size(in{1}(:,:,first_rand:end)));
            
            % adaptively shrink noise
            r = range(cat(3,in{:}),3); noisestd = .1 .* r;
            noisestd = repmat(noisestd,[1 1 size(in{1},3)]);
            
            
        elseif b == 1 && ~exist('genfun', 'var')  % then we do not have a custom input function, so do noise stuff
                    % add noise to inputs that are not best (1), if best = 1 (best of
        % last gen also, no change)
        % only works & is appropriate if param values are continuous (cont.
        % underlying space), which may not be the case with custom obj.
        % function
        
            switch paramtype
                case 'continuous'

                    in{1}(:,:,2:end) = in{1}(:,:,2:end) + noisestd(:,:,2:end) .* randn(size(in{1}(:,:,2:end)));

                    if exist('mininputs', 'var')
                        % we have a range of min/max values entered.  do not exceed max
                        % range
                        wh = in{1} < mininputs; in{1}(wh) = mininputs(wh);
                        wh = in{1} > maxinputs; in{1}(wh) = maxinputs(wh);

                    end

                case 'discrete'
                    % get new start vectors for 50%

                    wh = randperm(gensize); 
                    wh = wh(1 : round(gensize ./ 2));
                    wh(wh == 1) = [];
                    for j = wh
                        in{1}(:,:,j) = get_discrete_start_params(possiblevals, nvals, m, n);
                    end

                otherwise
                    error('Unknown parameter input type.');

            end

        end
        
        % --------------------------------------------------------------------
        % * print output, if requested
        % --------------------------------------------------------------------
        gentime(i) = etime(clock,t1);

        if doverbose
            erase_string(str);
            str = sprintf('Time: %3.0f ',etime(clock,t1)); fprintf(1,str);
        end

        % use original eff vector to get best, in case fitness is stochastic
        meff = eff(eff == max(eff)); meff = meff(1);
        beff(i) = meff;

        if doverbose
            fprintf(1,'Best: # %3.0f, Fitness: %3.2f',b,meff)
            figure(f1); plot(beff,'k','Linewidth',2); drawnow

            mybp = best_params{:};
            if size(mybp, 2) < length(mybp), mybp = mybp'; end
            all_bp(i, :) = mybp;
            
            figure(f2); cla; plot(all_bp); drawnow
            
            
            if length(mybp) < 10
                fprintf(1,'  Best params: ')
                fprintf(1,'%3.2f  ', mybp);        
            end

            if dostochastic && (b == 1 && ~exist('genfun', 'var'))  % changed 2010 to fix bug with no-stochastic option
                fprintf(' Noise added to param estimates with stdev %3.3f', noisestd(1));
            end

            fprintf('\n')

        end

        % --------------------------------------------------------------------
        % * Check for convergence
        % --------------------------------------------------------------------
        % two criteria: no efficiency above median efficiency, and no change in the last genconverge generations

        if i > genconverge
            last_gens = beff(i-genconverge+1:i);
            unique_in_last_gens = length(unique(last_gens));
            % if the above is one, there has been no change in last n
            % generations
        else
            unique_in_last_gens = Inf;
        end

        if ~(any(eff > median(eff))) || unique_in_last_gens == 1
            isconverged = 1;
            if doverbose, disp(['System converged at generation ' num2str(i)]), end

            beff(i+1:end) = [];
            gentime(i+1:end) = [];
            fit(:,i+1:end) = [];

            break
        end

    end  % End loop through generations

    if doverbose
        % final report
        ynstr = {'No' 'Yes'};
        fprintf(1,'\nGA Finished\n')
        fprintf(1,'___________________________________________________________\n')
        fprintf(1,'Generations: %3.0f \nOrganisms per generation: %3.0f\n',numgen,gensize);
        fprintf(1,'\nInitial fitness (average): %3.4f \nFinal fitness: %3.4f\n',mean(fit(:,1)),beff(end));
        fprintf(1,'\nConverged: %s \n',ynstr{isconverged+1});

        fprintf(1,'\nAverage time per iteration: %3.0f\n',mean(gentime));
        fprintf(1,'Average time per organism: %3.0f\n',mean(gentime)./gensize);
        fprintf(1,'Total time: %3.0f\n',etime(clock,t0));
        fprintf(1,'___________________________________________________________\n')
    end

    return





    % --------------------------------------------------------------------
    % * get random start state of gensize organisms
    % --------------------------------------------------------------------
    function in = random_start_state(inputs,gensize,paramtype)
        % determine range
        r = [min(inputs(:)) max(inputs(:))];

        [m, n] = size(inputs);

        if strcmp(paramtype, 'discrete')
            possiblevals = unique(inputs);  % should enter all possible values in input matrix
            nvals = length(possiblevals);
        end

        % create in variable = cells are inputs, columns are param sets, rows
        % params within sets.
        % first one is always the input you put in!
        in(:,:,1) = inputs; % + randn(size(inputs)) .* std(inputs(:));

        % rest of the population

        for j = 2:round(gensize./.8)

            switch paramtype
                case 'continuous'
                    % 80% is random within 2*range of inputs
                    tmp = rand(m,n); tmp = tmp.*r(2).*2 + r(1);

                case 'discrete'
                    tmp = get_discrete_start_params(possiblevals, nvals, m, n);
                otherwise
                    error('Unknown parameter distribution type');
            end

            in(:,:,j) = tmp;
        end

        % j = 3rd dim = organism

        for j = round(gensize./.8):gensize
            switch paramtype
                case 'continuous'
                    % 20% is input + noise
                    tmp = randn(size(inputs)) .* std(inputs(:)) + inputs;

                case 'discrete'
                    tmp = get_discrete_start_params(possiblevals, nvals, m, n);
            end

            in(:,:,j) = tmp;
        end

        return

    function tmp = get_discrete_start_params(possiblevals, nvals, m, n)
        tmp = zeros(m, n);
        for row = 1:m
            for col = 1:n
                % discrete value from possible set
                randp = randperm(nvals);
                tmp(row, col) = possiblevals(randp(1));
            end
        end

        return

    % --------------------------------------------------------------------
    % * get start state of gensize organisms that evenly spans param space
    % --------------------------------------------------------------------
    function [in,gensize] = sobol_start_state(inputs,gensize,doverbose)

        rows = size(inputs,1);
        cols = size(inputs,2);

        n = numel(inputs(:,:,1)); % n params, equals n dimensions

        npoints = ceil(gensize .^ (1./n)); % points in each dimension

        lowvals = inputs(:,:,1);
        hivals = inputs(:,:,2);

        str = 'combvals = combvec(';

        for i = 1:n
            % get start values for this parameter
            vals{i} = linspace(lowvals(i), hivals(i),npoints);

            if i > 1, str = [str ', ']; end
            str = [str 'vals{' num2str(i) '}'];
        end
        str = [str ');'];

        vals;
        % create combvec; each col is an organism, each row a
        % parameter
        combvals = [];
        eval(str)
        %%gensize = size(combvals,2); % new gen size; won't work for > 1
        %%inputs cell

        if doverbose
            fprintf(1,'Creating start state that spans the input range of param values.\n');
            fprintf(1,'Params: %3.0f, points in range for each param: %3.0f\n',n,npoints);
            fprintf(1,'nearest ideal generation size is %3.0f;  Actual gensize is %3.0f\n',size(combvals,2),gensize);
        end

        combvals = combvals(:,1:gensize);
                
        % reshape combvec to format of input

        for j = 1:gensize

            thisorg = reshape(combvals(:,j),rows,cols);
            in(:,:,j) = thisorg;

        end


        return


            

    % --------------------------------------------------------------------
    % * get best half, crossover to fill in lists
    % --------------------------------------------------------------------
function [newvec,b,best_params] = xover(paramvec,eff, paramtype)
    %
    % paramvec is 3-D matrix, 3rd dim is realization (organism)
    % columns are sets of params to be crossed over.
    %
    % b is index of best param vec - may have random noise, etc. added in
    % newvec
    % best_params is original best param vec, stored in newvec(:,:,1)

    % best one before xover
    b = find(eff == max(eff)); b = b(1); best_params = paramvec(:,:,b);


    w = find(eff > median(eff));

    % we can only do crossover if not all the designs are the same
    % ---------------------------------------------------
    if isempty(w) 
        
        warning('Extremely homogenous sample!')
        % add a little random noise to efficiency - 1% of var of efficiency
        eff = eff + randn(1,length(eff)) .* .01 * mean(eff);
        
        w = find(eff > median(eff));
        
        if isempty(w)
            w = 1:round(size(eff,2)./2);
        end

    end

    % save best half
    % ---------------------------------------------------
    newvec = paramvec(:,:,w);
    [nparams,dummy,last] = size(newvec);  % save size for adding random variation to existing


    %  %  CODE commented out for adding random noise right now.
    % % n = length(newvec(:));
    % % n = round(n/10);       % 10% - to add random noise

    % number of crossover points
    nxover = max(1,round(nparams ./ 50)); % 5;

    % fill in 2nd half with crossovers of first half
    % ---------------------------------------------------
    n_to_fill = size(paramvec,3) - last;

    for v = 1:n_to_fill

        % choose two random integers within 1st half
        w = ceil(rand(1,2) * last);

        babyv = [];

        for i = 1:size(newvec,2)        % for each set of parameters (columns)

            %babyv(:,i) = rcomb(newvec(:,i,w(1)),newvec(:,i,w(2)));

            babyv(:,i) = rcomb_multi(newvec(:,i,w(1)),newvec(:,i,w(2)),nparams,nxover);
        end

        % add to others
        newvec(:,:,last+v) = babyv;


    end

    % add random variation to existing best half
    % % % wh = round(rand(1,n) .* n);
    % % % wh(wh < 1) = 1; wh(wh > n) = n;
    % % %
    % % % nv = newvec(:,:,1:last);
    % % % nv(wh) = nv(wh) + randn(size(wh));
    % % % newvec = cat(3,nv,newvec(:,:,last+1:end));

    newvec(:,:,1) = best_params;      % re-insert best one


    return


    % --------------------------------------------------------------------
    % * crossover
    % --------------------------------------------------------------------
function c = rcomb(a,b)
    % combines 2 vectors of numbers at a random crossover point

    w = ceil(rand * (length(a) - 1));
    c = a(1:w);
    c = [c; b(w+1:end)];

    return

function  c = rcomb_multi(a,b,n,nxover)

    % % n = length(a);
    % % nxover = max(1,round(n ./ 50)); % 5;     % number of crossover points
    st = randperm(n);

    % divide the vectors up into chunks at random points
    st = [1 sort(st(1:nxover)) n];
    en = (st(2:end))-1; st = st(1:end-1); en(end) = n;

    c = a;
    for i = 1:2:length(st)
        c(st(i):en(i)) = b(st(i):en(i));
    end

    % example:
    % a = (1:100)';
    % b = 1000+a;
    % nxover = 10;
    % ...run code...
    % figure;plot(c)

    return


function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
    return



