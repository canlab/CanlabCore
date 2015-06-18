function varargout = igls_vs_weighted_vs_ols_sim1(meth, varargin)
    % igls_vs_weighted_vs_ols_sim1(meth, varargin)
    %
    % =====================================================================
    % Multi-function 'all in one' simulation tool for simulating group
    % results on regressions run within-subjects.
    %
    % A simple model, X, is fit to each of n_subjects data.
    % The predictors for this model are hard-coded as [1 0 1 0 1 0]'
    % The coefficients (betas) for each subject are subjected to a group
    % analysis treating subject as a random effect.  
    % The simulation compares the power and false positive rates for 
    % different analysis methods.
    % To accomplish this, null and alt hypothesis data is generated with
    % standard deviations and effect magnitudes specified by the user.
    % 
    % The simulation data is stored in a SIM object (structure), which is 
    % saved to disk and loaded when iterations are added or a summary is made.
    %
    % There are four main modes of operation for this function:
    % 'create'
    %       Creates a new SIM structure with sample sizes and variances
    %       specified by the user.  
    %       Critical values to enter are nobs_within, n_subjects,
    %       std_between, and std_within
    %       up to two of these can be vectors, and a full factorial
    %       combination of these levels will be created.
    %
    % 'addcell'
    %       Adds a single cell to the SIM object with samp. sizes and stds
    %       you specify.  Note: cell here does not refer to a cell in the
    %       Matlab sense, but rather a structure element in the SIM vector of
    %       structures that contains settings and results for one
    %       simulation with one combination of parameter settings.
    %
    %
    % 'iterations'
    %       Adds iterations to one or more of the cells in the SIM object.
    %       Iterations can be added by going to the directory where the SIM file
    %       is stored and entering the number to add when running the function.
    %
    % 'summary'
    %       Returns summary power and false positive rates from the SIM
    %       structure stored in the current directory, and creates maps of
    %       power and FPR by one or two outcome variables, depending on
    %       your specifications when you run it.
    % 
    % =====================================================================
    %
    % igls_vs_weighted_vs_ols_sim1('create', nobs_within, n_subjects, std_between, std_within)
    % igls_vs_weighted_vs_ols_sim1('create', [10 20 40 80 160], [5 10 20 40],1, 1);
    %
    % igls_vs_weighted_vs_ols_sim1('iterations', snum, iter)
    %
    % Tor Wager, Sept. 2007
    %
    % =====================================================================
    % 'summary' function
    % =====================================================================
    % [SUM, VALS] = igls_vs_weighted_vs_ols_sim1('summary', p-value threshold, SIM fields to organize results by);
    % [SUM, VALS] = igls_vs_weighted_vs_ols_sim1('summary', .05, {'n_subjects'});
    % Analysis names can be one name or two, or names followed by specific
    % levels.
    % Choices are: n_subjects, nobs_within, std_within, std_between
    %
    % [SUM, VALS] = igls_vs_weighted_vs_ols_sim1('summary', .05, {'std_within' 'std_between'});
    % saveas(gcf,'Summary_Power_FPR_plot','png');
    % [SUM, VALS] = igls_vs_weighted_vs_ols_sim1('summary', .05, {'std_within', 'std_between', [.1 .2]});
    %
    % Examples:
    % -------------------------------------------------------------------------
    % % Simulate effects of sample size within subjects and between subjects:
    % igls_vs_weighted_vs_ols_sim1('create', [10 20 40 80 160], [5 10 20 40], 1, 1);
    % igls_vs_weighted_vs_ols_sim1('iterations', 1, 100);
    % igls_vs_weighted_vs_ols_sim1('iterations', [], 1000);
    % [SUM, VALS] = igls_vs_weighted_vs_ols_sim1('summary', .05, {'std_within' 'std_between'});
    %

    switch meth

        case 'create'
            % Create a new SIM object and save in current directory
            % create_sim(nobs_within, n_subjects, std_between, std_within)
            create_sim(varargin{:})

        case 'addcell'
            % Add a cell to existing SIM object and re-save
            add_sim(varargin{:})

        case 'iterations'
            % Add iterations to existing SIM object
            % add_iterations(snum, iter)
            add_iterations(varargin{:})

        case 'summary'
            % Summarize results stored in SIM object
            [SUM, VALS] = summarize(varargin{:});

            varargout{1} = SUM;
            varargout{2} = VALS;

            %[SUM, VALS] = summarize(.05, {'n_subjects'})

    end


end   % End Main Function


% SUB-FUNCTIONS

function create_sim(nobs_within, n_subjects, std_between, std_within)
    % Enter row vectors

    SIM = [];
    if exist('SIM_igls_vs_weighted_vs_ols_sim1.mat', 'file')

        error('FILE: SIM_igls_vs_weighted_vs_ols_sim1.mat already exists!  Delete or go to a new directory to create.');
    else
        save SIM_igls_vs_weighted_vs_ols_sim1 SIM
    end

    mycell = {nobs_within, n_subjects, std_between, std_within};

    vals = combvec(mycell{:});

    nelements = size(vals, 2);

    diary('SIM_igls_vs_weighted_vs_ols_sim1.txt')

    fprintf('Creating new SIM structure with %3.0f SIM elements (conditions).\n They are:', nelements);

    print_matrix(vals, {}, {'nobs_within' 'nobs' 'std_between' 'std_within'})

    fprintf('\n');

    save SIM_igls_vs_weighted_vs_ols_sim1 -append SIM

    for i = 1:size(vals, 2)

        % add_sim(snum, nobs_within, n_subjects, std_between, std_within)
        add_sim(i, vals(1, i), vals(2, i), vals(3, i), vals(4, i))

    end

    diary off

end




function add_sim(snum, nobs_within, n_subjects, std_between, std_within)

    SIM = load_sim;

    % simulation element to add
    if isempty(snum)
        snum = length(SIM) + 1;
    end

    % Create structure
    % ------------------------------------------------------------------
    SIM(snum).null.true_intercept_mean = 0; % pop. mean (fixed effect)
    SIM(snum).null.true_slope_mean = 0;     % pop. mean slope (fixed effect)
    SIM(snum).null.true_intercept_std = 0;  % between-subjects variability in intercept
    SIM(snum).null.true_slope_std = 0;      % between-subjects variability in slope

    SIM(snum).nobs_within = nobs_within;
    SIM(snum).n_subjects = n_subjects;
    SIM(snum).std_within = std_within;

    SIM(snum).alt.true_intercept_mean = 1; % pop. mean (fixed effect)
    SIM(snum).alt.true_slope_mean = 1;     % pop. mean slope (fixed effect)
    SIM(snum).alt.true_intercept_std = std_between;  % between-subjects variability in intercept
    SIM(snum).alt.true_slope_std = std_between;      % between-subjects variability in slope

    SIM(snum).alt.effectsize_between = SIM(snum).alt.true_intercept_mean ./ std_between;
    SIM(snum).alt.effectsize_within = SIM(snum).alt.true_intercept_mean ./ std_within;

    if ~isfield(SIM(snum).null, 'unweighted_pvals'), SIM(snum).null.unweighted_pvals = []; end
    if ~isfield(SIM(snum).null, 'weighted_pvals'), SIM(snum).null.weighted_pvals = []; end
    if ~isfield(SIM(snum).null, 'igls_pvals'), SIM(snum).null.igls_pvals = []; end
    if ~isfield(SIM(snum).null, 'igls_p_randvariance'), SIM(snum).null.igls_p_randvariance = []; end

    if ~isfield(SIM(snum).alt, 'unweighted_pvals'), SIM(snum).alt.unweighted_pvals = []; end
    if ~isfield(SIM(snum).alt, 'weighted_pvals'), SIM(snum).alt.weighted_pvals = []; end
    if ~isfield(SIM(snum).alt, 'igls_pvals'), SIM(snum).alt.igls_pvals = []; end
    if ~isfield(SIM(snum).alt, 'igls_p_randvariance'), SIM(snum).alt.igls_p_randvariance = []; end

    SIM(snum).null.unweighted_time = [];
    SIM(snum).null.weighted_time = [];
    SIM(snum).null.igls_time = [];

    save SIM_igls_vs_weighted_vs_ols_sim1 -append SIM

end


function SIM = load_sim
    % ----------------------------------------------------------------
    % load SIM
    SIM = [];

    if exist('SIM_igls_vs_weighted_vs_ols_sim1.mat', 'file')
        load SIM_igls_vs_weighted_vs_ols_sim1
    else
        error('No File called SIM_igls_vs_weighted_vs_ols_sim1.mat found.  Create using ''create'' option.');
    end

end



function add_iterations(snum, iter)


    % load SIM
    % ----------------------------------------------------------------
    SIM = load_sim;


    % Batch mode: add to all SIM elements
    % ----------------------------------------------------------------
    if isempty(snum)
        for mysim = 1:length(SIM)
            add_iterations(mysim, iter);
        end

        return

    end
    % Continue on for 'single element' mode


    % last good iteration
    % ----------------------------------------------------------------
    % get starting iteration number

    % Find the last eligible good iteration; we'll add after this
    good = SIM(snum).null.unweighted_pvals ~= 0 & SIM(snum).null.weighted_pvals ~= 0 ...
        & SIM(snum).null.igls_pvals ~= 0 & SIM(snum).null.igls_p_randvariance ~= 0 ...
        & SIM(snum).alt.unweighted_pvals ~= 0 & SIM(snum).alt.weighted_pvals ~= 0 ...
        & SIM(snum).alt.igls_pvals ~= 0 & SIM(snum).alt.igls_p_randvariance ~= 0;

    if isempty(good)
        good = 0;
    else
        good = any(good);  % collapse across rows
        good = find(good);
    end

    if ~isempty(good), good = good(end); else good = 0; end

    fprintf(1,'In SIM(%3.0f) there are %3.0f valid iterations.  Adding %3.0f more.', snum, good, iter);

    startat = good + 1;

    % initialize iterations
    % ----------------------------------------------------------------
    for myfield = {'null' 'alt'};
        SIM(snum).(myfield{1}).unweighted_pvals(:, startat:end) = [];
        SIM(snum).(myfield{1}).unweighted_pvals = [SIM(snum).(myfield{1}).unweighted_pvals zeros(2, iter)];

        SIM(snum).(myfield{1}).weighted_pvals(:, startat:end) = [];
        SIM(snum).(myfield{1}).weighted_pvals = [SIM(snum).(myfield{1}).weighted_pvals zeros(2, iter)];

        SIM(snum).(myfield{1}).igls_pvals(:, startat:end) = [];
        SIM(snum).(myfield{1}).igls_pvals = [SIM(snum).(myfield{1}).igls_pvals zeros(2, iter)];

        SIM(snum).(myfield{1}).igls_p_randvariance(:, startat:end) = [];
        SIM(snum).(myfield{1}).igls_p_randvariance = [SIM(snum).(myfield{1}).igls_p_randvariance zeros(2, iter)];

    end

    SIM(snum).null.unweighted_time(:, startat:end) = [];
    SIM(snum).null.unweighted_time = [SIM(snum).null.unweighted_time  zeros(1, iter)];

    SIM(snum).null.weighted_time(:, startat:end) = [];
    SIM(snum).null.weighted_time = [SIM(snum).null.weighted_time  zeros(1, iter)];

    SIM(snum).null.igls_time(:, startat:end) = [];
    SIM(snum).null.igls_time = [SIM(snum).null.igls_time  zeros(1, iter)];






    % add iterations
    % ----------------------------------------------------------------
    for it = startat : startat + iter

        if mod(it - 1, 100) == 0 && it > startat
            fprintf(1, '%3.0f ', it);
            save SIM_igls_vs_weighted_vs_ols_sim1 -append SIM
        end

        x = zeros(SIM(snum).nobs_within, SIM(snum).n_subjects);
        x(1:2:end,:) = 1;          % create fixed x signal (predictor)

        % Generate Null data
        % ================================================
        myfield = 'null';

        % create individuals (with between-subjects noise) -- 2nd
        % level obs.
        c = normrnd(SIM(snum).(myfield).true_slope_mean, SIM(snum).(myfield).true_slope_std, SIM(snum).n_subjects, 1);              % slope between-n_subjectsjects variations
        d = normrnd(SIM(snum).(myfield).true_intercept_mean, SIM(snum).(myfield).true_intercept_std, SIM(snum).n_subjects, 1);      % intercept between-n_subjectsjects variations

        % Create y: Add between-n_subjects error (random effects) and measurement noise (within-n_subjectsjects error)
        for i = 1:SIM(snum).n_subjects
            y(:,i) = d(i) + c(i).*x(:,i) + normrnd(0, SIM(snum).std_within, SIM(snum).nobs_within, 1);
        end

        for i = 1:size(y, 2), YY{i} = y(:, i); end
        for i = 1:size(y, 2), XX{i} = x(:, i); end

        % Model fits, save p-values
        % ================================================
        t1 = clock;  stats = glmfit_multilevel(YY, XX, [], 'noverbose');
        SIM(snum).null.unweighted_time(it) = etime(clock, t1);
        SIM(snum).null.unweighted_pvals(:, it) = stats.p';

        % one-sample t-test, weighted by inv of btwn + within vars
        t1 = clock;  stats = glmfit_multilevel(YY, XX, [], 'noverbose', 'weighted');
        SIM(snum).null.weighted_time(it) = etime(clock, t1);
        SIM(snum).null.weighted_pvals(:, it) = stats.p';

        t1 = clock;  out = igls(y, x, 'noplot', 'noverbose');
        SIM(snum).null.igls_time(it) = etime(clock, t1);
        SIM(snum).null.igls_pvals(:, it) = out.p;
        SIM(snum).null.igls_p_randvariance(:, it) = out.p_randvariance;




        % Generate Alt data
        % ================================================
        myfield = 'alt';

        % create individuals (with between-subjects noise) -- 2nd
        % level obs.
        c = normrnd(SIM(snum).(myfield).true_slope_mean, SIM(snum).(myfield).true_slope_std, SIM(snum).n_subjects, 1);              % slope between-n_subjectsjects variations
        d = normrnd(SIM(snum).(myfield).true_intercept_mean, SIM(snum).(myfield).true_intercept_std, SIM(snum).n_subjects, 1);      % intercept between-n_subjectsjects variations

        % Create y: Add between-n_subjects error (random effects) and measurement noise (within-n_subjectsjects error)
        for i=1:SIM(snum).n_subjects
            y(:,i) = d(i) + c(i).*x(:,i) + normrnd(0, SIM(snum).std_within, SIM(snum).nobs_within, 1);
        end

        for i = 1:size(y, 2), YY{i} = y(:, i); end
        for i = 1:size(y, 2), XX{i} = x(:, i); end

        % Model fits, save p-values
        % ================================================
        stats = glmfit_multilevel(YY, XX, [], 'noverbose');
        SIM(snum).(myfield).unweighted_pvals(:, it) = stats.p';

        % one-sample t-test, weighted by inv of btwn + within vars
        stats = glmfit_multilevel(YY, XX, [], 'noverbose', 'weighted');
        SIM(snum).(myfield).weighted_pvals(:, it) = stats.p';

        out = igls(y, x, 'noplot', 'noverbose');
        SIM(snum).(myfield).igls_pvals(:, it) = out.p;
        SIM(snum).(myfield).igls_p_randvariance(:, it) = 2*out.p_randvariance;      % ML edit. This should be two-tailed.
    end

    fprintf(1, '\n')
    save SIM_igls_vs_weighted_vs_ols_sim1 -append SIM

end




function [SUM, VALS] = summarize(varargin)
    %[SUM, VALS] = summarize(.05, {'n_subjects'})

    pthresh = varargin{1};

    myplotvar = varargin{2};

    SIM = load_sim;

    %% Get FPR and TPR for each cell in simulation
    % --------------------------------------------------------------
    for mysim = 1:length(SIM)

        for myfield = {'null' 'alt'};

            % get names of fields
            names = fieldnames(SIM(mysim).(myfield{1})); whsave = false(size(names));
            for i = 1:length(names), if findstr(names{i}, '_p'), whsave(i) = 1; end, end
            names = names(whsave);

            for j = 1:length(names)
                % check iterations and remove invalid ones
                my_pvals = SIM(mysim).(myfield{1}).(names{j});

                invalid = any(my_pvals == 0, 1);

                if any(invalid)
                    fprintf('SIM(%03d).%s.%s Problem.\n', mysim, myfield{1}, names{j});
                    fprintf('Error: some p-values are exactly zero.  Invalid iterations? Removing %3.0f\n', sum(invalid));
                    my_pvals(:, invalid) = [];
                end

                % get p-values
                switch myfield{1}
                    case 'null'
                        field_save_name = [names{j} '_fpr'];
                    case 'alt'
                        field_save_name = [names{j} '_pow'];
                    otherwise
                        error('Uh-oh! Illegal field name.')
                end

                SIM(mysim).(myfield{1}).(field_save_name) = sum(my_pvals < pthresh, 2) ./ size(my_pvals, 2);

            end

        end

    end

    % Create summary object with overall results
    % --------------------------------------------------------------
    SUM.nobs_within = cat(2, SIM(:).nobs_within);
    SUM.n_subjects = cat(2, SIM(:).n_subjects);
    SUM.std_within = cat(2, SIM(:).std_within);

    for i = 1:length(SIM), SUM.std_between(1, i) = SIM(i).alt.true_intercept_std;  end

    for i = 1:length(SIM)
        SUM.unweighted_fpr(:,i) = SIM(i).null.unweighted_pvals_fpr;
        SUM.weighted_fpr(:,i) = SIM(i).null.weighted_pvals_fpr;
        SUM.igls_fpr(:,i) = SIM(i).null.igls_pvals_fpr;
        SUM.igls_p_randvariance_fpr(:,i) = SIM(i).null.igls_p_randvariance_fpr;

        SUM.unweighted_pow(:,i) = SIM(i).alt.unweighted_pvals_pow;
        SUM.weighted_pow(:,i) = SIM(i).alt.weighted_pvals_pow;
        SUM.igls_pow(:,i) = SIM(i).alt.igls_pvals_pow;
        SUM.igls_p_randvariance_pow(:,i) = SIM(i).alt.igls_p_randvariance_pow;

    end

    % GET averages by levels of selected variables
    % --------------------------------------------------------------
    %myplotvar = input cell array     {'nobs_within'};

    % get levels of vars to plot
    VALS.myvars = {};
    VALS.levels = {};

    for i = 1:length(myplotvar)
        if ischar(myplotvar{i})

            VALS.myvars{end+1} = myplotvar{i};

            if length(myplotvar) > i && ~ischar(myplotvar{i + 1})
                VALS.levels{end+1} = myplotvar{i + 1};
            else
                VALS.levels{end+1} = unique(SUM.(myplotvar{i}));
            end

            VALS.(myplotvar{i}) = VALS.levels{end};
        end
    end

    switch length(VALS.myvars)
        case 1
            % one-dim plot, by one var

            for j = 1:length(VALS.levels{1})

                wh = SUM.(VALS.myvars{1}) == VALS.levels{1}(j);  % which SIM indices are of this level of chosen variable

                for name = {'weighted_fpr' 'unweighted_fpr' 'igls_fpr' 'igls_p_randvariance_fpr'};

                    VALS.(name{1})(j) = mean(SUM.(name{1})(2, wh));

                end

                for name = {'weighted_pow' 'unweighted_pow' 'igls_pow' 'igls_p_randvariance_pow'};

                    VALS.(name{1})(j) = mean(SUM.(name{1})(2, wh));

                end

            end

        case 2
            % two-dim plot, by two vars

            for j = 1:length(VALS.levels{1})

                for k = 1:length(VALS.levels{2})

                    % which SIM indices are of this level of chosen variables
                    wh = SUM.(VALS.myvars{1}) == VALS.levels{1}(j) & SUM.(VALS.myvars{2}) == VALS.levels{2}(k);

                    for name = {'weighted_fpr' 'unweighted_fpr' 'igls_fpr' 'igls_p_randvariance_fpr'};

                        VALS.(name{1})(j, k) = mean(SUM.(name{1})(2, wh));

                    end

                    for name = {'weighted_pow' 'unweighted_pow' 'igls_pow' 'igls_p_randvariance_pow'};

                        VALS.(name{1})(j, k) = mean(SUM.(name{1})(2, wh));

                    end

                end

            end



        otherwise
            error('Enter at most 2 names in vars to plot');
    end

    % Create plot
    % --------------------------------------------------------------

    colors = {'ro-' 'gs-' 'bv-' 'cd-' 'm^-' 'kx-' 'y+-'};

    switch length(VALS.myvars)

        case 1

            create_figure('Simulation Results', 1, 2);

            cindx = 1;
            subplot(1, 2, 1);

            names = {'weighted_fpr' 'unweighted_fpr' 'igls_fpr' 'igls_p_randvariance_fpr'};
            for name = names
                plot(VALS.(myvars{1}), VALS.(name{1}), colors{cindx}, 'LineWidth', 2, 'MarkerFaceColor', colors{cindx}(1));
                cindx = cindx + 1;
            end

            title('False Positive Rate');
            xlabel(VALS.myvars{1});


            cindx = 1;
            subplot(1, 2, 2);

            names = {'weighted_pow' 'unweighted_pow' 'igls_pow' 'igls_p_randvariance_pow'};
            for name = names
                plot(VALS.(VALS.myvars{1}), VALS.(name{1}), colors{cindx}, 'LineWidth', 2, 'MarkerFaceColor', colors{cindx}(1));
                cindx = cindx + 1;
            end

            title('Power');
            xlabel(VALS.myvars{1});

            legend(names);

        case 2

            create_figure('Simulation Results', 2, 4);

            analysisname = {'unweighted_fpr' 'weighted_fpr' 'igls_fpr' 'igls_p_randvariance_fpr' 'unweighted_pow' 'weighted_pow' 'igls_pow' 'igls_p_randvariance_pow'};
            for an = 1:length(analysisname)

                subplot(2, 4, an);
                for i = 1:length(VALS.levels{2})
                    plot(VALS.levels{1}, VALS.(analysisname{an})(:, i), colors{i}, 'LineWidth', 2, 'MarkerFaceColor', colors{i}(1));
                end

                xlab = VALS.myvars{1};
                if ~isempty(findstr(analysisname{an}, 'fpr'))
                    ylab = 'False Positive Rate';
                    ylimit = [0 min(max(pthresh*4, .05), 1)];

                    str = analysisname{an};  str(str == '[' | str == ']') = []; str(str == '_') = ' '; str = str(1:end-4);
                    title(str, 'FontSize', 24);

                    plot_horizontal_line(pthresh);

                elseif ~isempty(findstr(analysisname{an}, 'pow'))
                    ylab = 'Detection Power';
                    ylimit = [0 1];
                end

                axis_labels(xlab, ylab);
                set(gca,'YLim', ylimit);

                % legend
                if an == 1
                    str = mat2str(VALS.levels{2}); str(str == '[' | str == ']') = []; str = parse_char_to_cell(str, 'space');
                    legend(str)
                end

            end

    end   % one-d vs two-d plot

    scn_export_papersetup(700);

end



function axis_labels(xlab, ylab)

    xlab(1) = upper(xlab(1)); xlab(xlab == '_') = ' ';
    xlabel(xlab);

    ylab(1) = upper(ylab(1)); ylab(ylab == '_') = ' ';
    ylabel(ylab);

end

