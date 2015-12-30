function varargout = scnlab_outlier_id(varargin)
% Methods (modes of operation)
%
% :Setup:
%
%      Run this method first to generate an options structure OPT
%      that can be passed in along with any data vector for speedy
%      processing
%
% :Usage:
% ::
%
%     OPT = scnlab_outlier_id('setup', 'tr', 2, 'spersess',
%                             [184 184 184 184 184 184], 'dummy', 1:3,
%                             'hp', 100, 'mad', 4, 'niter', 3, 'mvmt', mvmt);
%
% :Data:
%      Run this method second with an already-created OPT
%      and a data vector from one time series
% ::
%
%     [y2, outliers, num_outliers, mvmt_rsquare] = scnlab_outlier_id('data', y, 'options', OPT);
%
% all outputs:
% ::
%
%     [y2, out, nout, mvmtrsq, mvmt_baseline_rsquare, yperc, rawvarp, rawvarF, percvarp, percvarF, ...
%     ybase] = scnlab_outlier_id('data', y, 'options', OPT);
%
% :Example:
% ::
%
%    [dat, volInfo] = iimg_get_data('graymask.img', imgs);
%    y = dat(:,1);
%    % SETUP:
%    OPT = scnlab_outlier_id('setup', 'tr', 2, 'spersess', [184 184 184 184 184 184], 'dummy', 1:3, 'hp', 100, 'mad', 4, 'niter', 3, 'mvmt', mvmt);
%    % RUN:
%    [y2, outliers, num_outliers] = scnlab_outlier_id('data', y, 'options', OPT);
%
%    % run on whole brain
%    [dat, volInfo] = iimg_get_data('graymask.img', imgs);
%    OPT = scnlab_outlier_id('setup', 'tr', 2, 'spersess', [184 184 184 184 184 184], 'dummy', 1:2, 'hp', 100, 'mad', 4, 'niter', 5, 'mvmt', mvmt);
%    OPT.doplot = 0;
%    OPT.verbose = 0;
%    fhandle = @(y) scnlab_outlier_id('data', y, 'options', OPT);
%    y2 = matrix_eval_function(dat, fhandle)';
%
% ..
%    Created June 2007
%    Tor Wager
% ..


    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'setup', meth = 'setup';
                case 'data', meth = 'data';
            end
        end
    end

    OPT.tr = [];
    OPT.spersess = [];
    OPT.dummy = [];
    OPT.hp = [];
    OPT.mad = 4;
    OPT.niter = 3;
    OPT.doplot = 1;
    OPT.verbose = 1;
    OPT.mvmt = [];
    OPT.dopercent = 1;

    switch meth

        % ---------------------------------------------------------------------
        % *
        % *
        % *  SETUP
        % *     Run this method first to generate an options structure OPT
        %       that can be passed in along with any data vector for speedy
        %       processing
        % *
        % ---------------------------------------------------------------------

        case 'setup'


            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}

                        % functional commands
                        case 'tr', OPT.tr = varargin{i+1};
                        case 'spersess', OPT.spersess = varargin{i+1};
                        case 'dummy', OPT.dummy = varargin{i+1};
                        case 'hp', OPT.hp = varargin{i+1};
                        case 'mad', OPT.mad = varargin{i+1};
                        case 'niter', OPT.niter = varargin{i+1};
                        case 'plot', OPT.doplot = varargin{i+1};
                        case 'verbose', OPT.verbose = varargin{i+1};

                        case 'mvmt', OPT.mvmt = varargin{i+1};

                        case 'nopercent', OPT.dopercent = 0;

                        case {'setup', 'data'} % do nothing

                        otherwise, warning(['Unknown input string option:' varargin{i}]);
                    end
                end
            end


            if OPT.verbose
                disp('Checking inputs.')
            end

            if isempty(OPT.tr), error('Enter ''tr'' followed by tr of acquisition.'); end
            if isempty(OPT.spersess), error('Enter ''spersess'' followed by vector of number of images in each session.'); end

            OPT.spersess = OPT.spersess(:)'; %enforce row vector
            
            % set up matrices for linear, invariant part
            % ---------------------------------------------------------------------
            nimgs = sum(OPT.spersess);
            nsess = length(OPT.spersess);

            y = ones(nimgs, 1);
            OPT.descrip = 'Fields below will be applied to data during processing.';
            OPT.descrip2 = 'They contain high-pass filter, session mean, and dummy scan information';
            % ***** does most of the work ******
            [y, OPT.IpinvI, OPT.HPmatrix] = hpfilter(y, OPT.tr, OPT.hp, OPT.spersess,[],OPT.dummy);


            % Get problematic images from motion parameters
            % ---------------------------------------------------------------------
            if ~isempty(OPT.mvmt)

                nparams = size(OPT.mvmt, 2);

                if size(OPT.mvmt, 1) ~= length(y), warning('Length of movement parameters not equal to number of images!!!'); end

                OPT.mvmt_diff = [zeros(1, nparams); diff(OPT.mvmt)];

                OPT.moved_more_than_point1 = any(abs(OPT.mvmt_diff(:,4:6)) > .1, 2);
                OPT.rotate_more_than_point002 = any(abs(OPT.mvmt_diff(:,1:3)) > .002, 2);

                % convert to MADs
                OPT.mvmt_diff = OPT.mvmt_diff ./ repmat(mad(OPT.mvmt_diff), size(OPT.mvmt_diff,1), 1);

                % ID outliers/images with high movement

                OPT.mvmt_outliers = any(abs(OPT.mvmt_diff) > OPT.mad, 2);

                % these would mark vals on either side of rapid transition
                %OPT.mvmt_outliers + [zeros(1, nparams); OPT.mvmt_outliers(1:end-1, :)];

                % Outliers: have to have substantial motion, and also be
                % out of subject's distribution based on MAD
                OPT.mvmt_outliers = any(OPT.mvmt_outliers, 2) & (OPT.moved_more_than_point1 | OPT.rotate_more_than_point002);

                % other movement things of interest
                OPT.within_sess_mvmt = OPT.mvmt - OPT.IpinvI * OPT.mvmt;

                OPT.total_mvmt_displacement = max(OPT.mvmt) - min(OPT.mvmt);
                en = cumsum(OPT.spersess);
                st = [1 en(1:end-1)+1];

                for jj = 1:length(OPT.spersess)
                    wh = st(jj):en(jj);
                    OPT.within_mvmt_displacement(jj, :) = max(OPT.within_sess_mvmt(wh,:)) - min(OPT.within_sess_mvmt(wh,:));
                end

                if OPT.doplot
                    create_figure1
                    nfigrows; nfigcols;
                    plot_movement
                end

            end


            % print output

            if OPT.verbose
                fprintf(1,'\nscnlab_outlier_id setup\n----------------------------------------\n');
                fprintf('Data: %3.0f observations in %3.0f runs (sessions) \n', nimgs, nsess);
                fprintf('\t Images in each run:')
                fprintf('%3.0f  ', OPT.spersess);
                fprintf('\n');
                fprintf('TR: %3.2f\n', OPT.tr);

                nystr = {'No' 'Yes'};

                fprintf('High-pass filter: %s ', nystr{(~isempty(OPT.hp)) + 1});
                if ~isempty(OPT.hp)
                    fprintf(', HP length: %3.2f\n', OPT.hp);
                else
                    fprintf('\n');
                end


                fprintf('Movement parameters entered: %s \n', nystr{(~isempty(OPT.mvmt)) + 1});
                %                 if ~isempty(OPT.mvmt)
                %                     fprintf('HP length: %3.2f\n', OPT.hp);
                %                 else
                %                     fprintf('\n');
                %                 end

                fprintf('MAD threshold for outliers: %3.2f\n', OPT.mad);
                fprintf('Number of iterations: %3.0f\n', OPT.niter);

                if ~isempty(OPT.mvmt)
                    fprintf(1,'\nMovement summary\n----------------------------------------\n');
                    fprintf('Maximum displacement across experiment\n')
                    fprintf('r\tp\ty\tx\ty\tz\t\n')
                    fprintf('%3.2f\t', OPT.total_mvmt_displacement);
                    fprintf('\n\n');

                    fprintf('Maximum displacement within sessions\n')
                    fprintf('r\tp\ty\tx\ty\tz\n')
                    for jj = 1:length(OPT.spersess)
                        fprintf('%3.2f\t', OPT.within_mvmt_displacement(jj, :));
                        fprintf('\n');
                    end
                    fprintf('(mean within session displacement)\n')
                    fprintf('%3.2f\t', mean(OPT.within_mvmt_displacement));
                    fprintf('\n\n');

                    fprintf('Images with displacement > 0.1 mm: %3.0f\n', sum(OPT.moved_more_than_point1))
                    fprintf('Images with rotation > 0.002 radians: %3.0f\n', sum(OPT.rotate_more_than_point002))
                    fprintf('Images excluded as motion transients: %3.0f\n', sum(OPT.mvmt_outliers))

                end

            end

            varargout{1} = OPT;













            % ---------------------------------------------------------------------
            % *
            % *
            % *  DATA
            % *     Run this method second with an already-created OPT
            %       and a data vector from one time series
            %
            % *
            % ---------------------------------------------------------------------

        case {'data'}
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}

                        % functional commands
                        case {'opt', 'options'}, OPT = varargin{i+1};
                        case {'y', 'data'}, y = varargin{i+1};

                        case 'plot', OPT.doplot = varargin{i+1};
                        case 'verbose', OPT.verbose = varargin{i+1};

                        case {'setup'} % do nothing

                        otherwise, warning(['Unknown input string option:' varargin{i}]);
                    end
                end
            end

            n = size(y,1);
            outliers = false(size(y));

            % add movement outliers, if we have them
            % ---------------------------------------------------------------------
            if ~isempty(OPT.mvmt)
                outliers(OPT.mvmt_outliers) = 1;
            end


            if OPT.doplot
                % ---------------------------------------------------------------------
                create_figure1
                subplot(nfigrows, nfigcols, 1);

                if ~isempty(OPT.mvmt)
                    plot_movement
                end

            end


            % setup for moving average
            % starting and ending images for each session
            st = cumsum([1 OPT.spersess(1:end-1)]);
            en = cumsum(OPT.spersess);
            yfit_baseline = zeros(n, 1);

            for i = 1:OPT.niter

                % moving average method
                % ---------------------------------------------------------------------
                iteration_moving_average

                % iteration_hp_filter

            end


            % final step
            % interpolate outliers
            % ---------------------------------------------------------------------

            y = interpolate_outliers(y, yfit_baseline, outliers, n);
            


            yperc = 0;
            raw_equalvar_p = 0;
            raw_equalvar_F = 0;
            perc_equalvar_p = 0;
            perc_equalvar_F = 0;

            if OPT.dopercent
                yperc = 100 .* y ./ yfit_baseline;

                % test equality of variance for y and yperc
                % if percent change is better scaling, variances should be
                % more equal across groups

                for ss = 1:length(st)
                    
                    group(st(ss):en(ss)) = ss;
                    
                end

                % test equality of variances across sessions
                [raw_equalvar_p, stats] = vartestn(y, group, 'off', 'robust');
                raw_equalvar_F = stats.fstat;

                [perc_equalvar_p, stats] = vartestn(yperc, group, 'off', 'robust');
                perc_equalvar_F = stats.fstat;

            end
            
            varargout{6} = yperc;
            varargout{7} = raw_equalvar_p;
            varargout{8} = raw_equalvar_F;
            varargout{9} = perc_equalvar_p;
            varargout{10} = perc_equalvar_F;


            if OPT.doplot
                plot_this_iteration

                subplot(nfigrows, nfigcols, 7)

                %  plot sessions
                draw_session_boxes(OPT, y)
            end


            % movement-related stats
            % ---------------------------------------------------------------------
            mvmt_rsquare = 0; vy = var(y);
            
            if ~isempty(OPT.mvmt) && ~isnan(vy) && vy > 0
                
                [b, bint, r, rint, stats] = regress(y, [OPT.mvmt ones(n, 1)]);
                mvmt_rsquare = stats(1);
                
                [b, bint, r, rint, stats] = regress(yfit_baseline, [OPT.mvmt ones(n, 1)]);
                mvmt_baseline_rsquare = stats(1);
                
            end


            varargout{1} = y;
            varargout{2} = outliers;
            varargout{3} = sum(outliers);
            varargout{4} = mvmt_rsquare;
            varargout{5} = mvmt_baseline_rsquare;
            varargout{11} = yfit_baseline;          % at end, so we don't have to save lots of data in matrix_eval_fcn


            if OPT.verbose
                fprintf(1,'\nData summary\n----------------------------------------\n');
                fprintf('Outliers\n')
                fprintf('Total: %3.0f\n', sum(outliers));
                fprintf('Data for outliers estimated using linear interpolation.\n')
                if ~isempty(OPT.mvmt)
                    fprintf('Movement-related exclusions: %3.0f\n', sum(OPT.mvmt_outliers));
                    fprintf('% Var explained by mvmt params: %3.2f\n', mvmt_rsquare);
                end

                fprintf('\n');

            end

        otherwise
            error('Unknown method')

    end




    % ---------------------------------------------------------------------
    % INLINE FUNCTIONS
    % ---------------------------------------------------------------------

    function create_figure1

        nfigrows = 4;
        nfigcols = 2;

        create_figure('Data Detail', nfigrows, nfigcols);

        subplot(nfigrows,nfigcols,3);
        %  plot sessions
        draw_session_boxes(OPT, y)

        plot(y, 'k', 'LineWidth', 2);
        title('Original data')
        set(gca,'XLim',[0 length(y)])


        subplot(nfigrows,nfigcols,4);
        shaded_hist(y);
        title('Histogram');

        drawnow

    end

    function plot_movement

        subplot(nfigrows,nfigcols,1);

        %  plot sessions
        draw_session_boxes(OPT, OPT.mvmt)
        plot(OPT.mvmt)

        title('Movement parameters')
        ylabel('mm or degrees');
        set(gca,'XLim',[0 length(y)])


        subplot(nfigrows,nfigcols,2);
        plot(OPT.mvmt_diff)

        %  plot sessions
        draw_session_boxes(OPT, OPT.mvmt)

        wh_outliers = find(OPT.mvmt_outliers);
        if(~isempty(wh_outliers))
            for jj = wh_outliers
                han = plot_vertical_line(jj,'k'); 
                set(han,'Color',[.7 .7 .7]); 
            end
        end
        plot(OPT.mvmt_diff);
        %plot(repmat(find(OPT.mvmt_outliers), 1, size(OPT.mvmt, 2)), OPT.mvmt_diff(OPT.mvmt_outliers, :), 'rs');

        title('Displacement (MADs)')
        ylabel('Med Abs Dev');
        set(gca,'XLim',[0 length(y)], 'YLim', [-OPT.mad * 4 OPT.mad * 4])

    end



    function iteration_moving_average



        % remove drift and session effects
        % ---------------------------------------------------------------------
        % for each session
        for ss = 1:length(st)

            yfit_baseline(st(ss):en(ss)) = moving_average('gaussian',y(st(ss):en(ss)), round(OPT.hp ./ OPT.tr));

        end

        % get adjusted timeseries
        yadj = y - yfit_baseline;

        if OPT.doplot

            if i == 1, plot_first_pass; end
            plot_this_iteration

        end

        % ID outliers on residuals and update outlier list
        % ---------------------------------------------------------------------
        mady = mad(yadj);
        wh = abs(yadj) > mady * OPT.mad;
        outliers = outliers | wh;

        % interpolate at outlier values
        y(outliers) = yfit_baseline(outliers);

    end




    function iteration_hp_filter
        % impute mean to identified outliers
        % ---------------------------------------------------------------------
        y(outliers) = repmat(mean(y, 1), sum(outliers), 1);  % impute mean

        % remove drift and session effects
        % ---------------------------------------------------------------------
        y = hpfilter(y, [], OPT.HPmatrix, OPT.spersess, OPT.IpinvI);

        if OPT.doplot

            if i == 1, plot_first_pass; end
            plot_this_iteration

        end

        % ID outliers on residuals and update outlier list
        % ---------------------------------------------------------------------
        mady = mad(y); % repmat(mad(y), n, 1);
        wh = (abs(y) > OPT.mad * mady);

        outliers = outliers | wh;                            % add to list of outliers
    end


    function plot_this_iteration
        subplot(nfigrows, nfigcols, 7)
        cla
        if exist('yadj', 'var')     % moving average
            myploty = yadj;
        else   % hpfilter method
            myploty = y;
        end

        plot(myploty, 'k');
        hold on;
        plot(find(outliers), myploty(outliers), 'rs', 'LineWidth', 2);

        set(gca,'XLim',[0 length(myploty)])
        title(['Iteration ' num2str(i)]);

        subplot(nfigrows,nfigcols,8);
        cla
        shaded_hist(myploty);
        title('Histogram');

        drawnow
        pause(.5)
    end


    function plot_first_pass
        subplot(nfigrows, nfigcols, 5)
        cla
        plot(y, 'k');
        hold on;
        plot(find(outliers), y(outliers), 'rs', 'LineWidth', 2);

        set(gca,'XLim',[0 length(y)])
        title('First pass');

        if exist('yfit_baseline', 'var') % moving average only
            plot(yfit_baseline, 'r');
        end

        subplot(nfigrows,nfigcols,6);
        cla
        shaded_hist(y);
        title('Histogram');

        drawnow
        pause(.5)
    end

end         % END MAIN FUNCTION



function draw_session_boxes(OPT, y)


    nsess = length(OPT.spersess);
    colors = {'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k'};  %{'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
    while nsess > length(colors), colors = [colors colors]; end

    yval = min(y(:)) + .05 * min(y(:));
    yheight = (max(y(:)) + .05 * max(y(:))) - min(y(:));
    c = cumsum(OPT.spersess);
    st = [0 c(1:end-1)];

    for jj = 1:2:nsess
        h1 = drawbox(st(jj),OPT.spersess(jj),yval,yheight,colors{jj}(1));
        set(h1, 'FaceAlpha', .10, 'EdgeColor', 'none');
    end

end


% -------------------------------------------------------------------------
% Create a shaded gray-scale histogram with .05 2-tailed in darker gray
% -------------------------------------------------------------------------
function shaded_hist(a, xx)
    if nargin < 2
        nbins = max(10, round(length(a) ./ 20));
        nbins = min(nbins, 1000);
        [h, xx] = hist(a, nbins);
    else
        h = hist(a, xx);
    end
    han = bar(xx, h);
    set(han, 'FaceColor', [.7 .7 .7]);
    set(han, 'EdgeColor', [.7 .7 .7]);
    wh = (xx > prctile(a, 97.5) | xx < prctile(a, 2.5));
    h(~wh) = 0;
    if any(wh')
        hold on;
        han = bar(xx, h);
        set(han, 'FaceColor', 'r', 'EdgeColor', [.3 .3 .3]);
    end
    plot_vertical_line(0);
end


function y = interpolate_outliers(y, yfit_baseline, outliers, n)
    y = y - yfit_baseline;

    time = 1:n;
    whout = find(outliers);


    yinterp = interp1(time(~outliers),y(~outliers),whout, 'linear');

        % NaNs are used at ends
    mynans = isnan(yinterp);
    yinterp(mynans) = 0;   % yfit_baseline(whout(mynans));

    
    y(outliers) = yinterp;
    % %
    % %             % remove drift and session effects
    % %
    % %             y = hpfilter(y, [], OPT.HPmatrix, OPT.spersess,
    % OPT.IpinvI);
end
            
