% function varargout = trial_level_beta3([onsets or design matrix], [data vector])
%
% This function can do three things:
% 1) Build a trial-level design matrix X from a vector of onsets
% 2) Fit a design matrix X to data to get trial-level betas
% 3) Estimate height, time to peak, and width across individual trials
%
% Tor Wager
% Last modified: 4/1/13 (Luka)
%
% output types are:
% 'betas', 'design'
%
% Options:
%   'tr', tr_in_seconds (DEFAULT: 2)
%   'stimlength', stimulus_length_in_seconds (DEFAULT: 10) (only matters to painX)
%
% Basis sets:
% 'basistype' followed by basis set name
% choices:
% 'spm+disp' : spm + time and disp derivatives
% 'spmgamma2': two gamma functions, second delayed (10 s peak)
% 'hrf' : spm canonical HRF
% 'gamma3' : three SPM HRFs, the original + shifted forwards and back one TR
% 'gamma4' : four canonical HRFs, orig, back one TR, forward 1 and 2 TRs
% 'pain10' : 3-param BF designed for TR=2, NSF and expectancy pain studies (stimlength 10s)
% 'painX' : 3-param BF constructed similar to pain10 but with variable TR
%             and stimlength
% 'flat10' :  3-param BF designed for TR=2, des. for NSF and expectancy anticipation
%
% Examples:
% Data fitting and model building:
% [betas, f, X, px, bf] = trial_level_beta3('onsets', onsets, 'rows', 6880, 'output', 'betas', 'plot', 1);
%
% Data fitting, given model as input
% [betas, f, X2, px2] = trial_level_beta3('pinvx', px, 'X', X, 'output', 'betas', 'plot', 1, 'data', y);
% [betas, f, X2, px2] = trial_level_beta3('X', X, 'output', 'betas', 'plot', 1, 'data', y);
%
% Data fitting, providing trial-level height, time-to-peak, and width as
% output:
% [betas, f, X2, px2, bf, h, t, w] = trial_level_beta3('X', X, 'pinvx', px, 'bf', bf, 'output', 'betas', 'plot', 1, 'data', y);
%
% Design building
% [X, px] = trial_level_beta3('onsets', onsets, 'output', 'design');
% [X, px, bf] = trial_level_beta3('onsets', onsets, 'rows', rowsz, 'output', 'design', 'basistype', 'gamma4');
%
% Only get basis set:
% [X, px, bf] = trial_level_beta3('onsets', [], 'output', 'design', 'basistype', 'pain10');
%
% Example: pain: build with plot, then fit with model as input (faster)
% [betas, f, X, px, bf] = trial_level_beta3('data', mm, 'onsets', ons, 'rows', 640, 'output', 'betas', 'plot', 1);
% [betas, f, X2, px2, bf, h, t, w] = trial_level_beta3('X', X, 'pinvx', px, 'bf', bf, 'output', 'betas', 'plot', 0, 'data', mm);
%
% Build a design with TR=2 and a different basis set, then estimate:
% [trialX, trial_px, bf] = trial_level_beta3('onsets', times, 'output', 'design', 'rows', 162, 'basistype', 'gamma3', 'tr', 2);
% [betas, f, X2, px2, bf, h, t, w, auc] = trial_level_beta3('X', trialX, 'pinvx', trial_px, 'bf', bf, 'output', 'betas', 'plot', 1, 'data', datv);
%
% Specify a certain number of trials for h, t, w (i.e., when adding nuisance params to trialX)
% [betas, f, X2, px2, bf, h, t, w, auc] = trial_level_beta3('X', trialX, 'pinvx', trialpx, 'bf', bf, 'output', 'betas', 'plot', 1, 'data', datv, 'ntrials', 20);

function varargout = trial_level_beta3(varargin)

    rowsz = [];
    y = [];
    onsets = [];
    px = [];
    X = [];
    htw = 1;
    doplot = 0;
    output = 'betas';
    tr = 2;
    stimlength = 10;
    basistype = 'spm+disp';
   

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'betas'
                case 'design'

                    % functional commands
                case 'rows', rowsz = varargin{i+1};
                case 'data', y = varargin{i+1};
                case 'onsets', onsets = varargin{i+1};
                case 'pinvx', px = varargin{i+1};
                case 'X', X = varargin{i+1};
                case 'bf', bf = varargin{i+1};
                case 'plot', doplot = varargin{i+1};
                case 'output', outtype = varargin{i+1};
                case 'htw', htw = varargin{i+1};
                case 'ntrials', ntrials = varargin{i+1};

                    % parameters
                case 'stimlength', stimlength = varargin{i+1};
                case 'tr', tr = varargin{i+1};
                case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end


    switch outtype
        case 'betas'
            % ---------------------------------------------------------------------
            % Create beta series by fitting model, and plot if requested
            % ---------------------------------------------------------------------
            if isempty(y), error('You need data to get betas output. Enter ''data'', y'); end

            rowsz = length(y);

            if isempty(px)
                if ~isempty(X), px = pinv(X);
                else
                    % we need to create X from onsets
                    [X, px, bf] = trial_level_beta3('onsets', onsets, 'rows', rowsz, 'output', 'design', 'plot', doplot);
                end
            end

            b = px * y;
            f = X * b;

            varargout{1} = b; varargout{2} = f;
            if nargout > 1, varargout{2} = f; end
            if nargout > 2, varargout{3} = X; end
            if nargout > 3, varargout{4} = px; end
            if nargout > 4, varargout{5} = bf; end

            if htw || nargout > 4
                % get height, time to peak, and width for each
                if exist('ntrials', 'var')
                    [h, t, w, auc] = get_beta_htw(b(1:(ntrials*size(bf, 2))), doplot, bf);
                else
                    [h, t, w, auc] = get_beta_htw(b, doplot, bf);
                end

                varargout{6} = h;
                varargout{7} = t;
                varargout{8} = w;
                varargout{9} = auc;

                % TRIAL AMPLITUDE PLOT
                if doplot
                    create_figure('Trial beta timeseries', 2, 1);
                    subplot(2, 1, 2);
                    
                    plot(h, 'ko', 'MarkerFaceColor', [.5 .5 .5]);
                    
                    plot_horizontal_line(0, 'k');
                    xlabel('Trial');
                    ylabel('Amplitude');
                    title('Trial-level amplitude estimates');
                end
            end


            if doplot
                create_figure('Trial beta timeseries', 2, 1, 1);
                subplot(2, 1, 1);
                
                plot(y, 'k');
                plot(f, 'r');
                legend({'Data' 'Fitted'});

                if exist('onsets', 'var') && ~isempty(onsets)
                    plot(onsets+t, h+b(end), 'ko', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [.8 .8 .8]);
                    hh = plot_onsets(onsets, 'k', -10, 5, 1);  % last arg is length in TRs
                end
                drawnow;
            end


        case 'design'
            % ---------------------------------------------------------------------
            % Build design matrix from onsets
            % ---------------------------------------------------------------------

            % CURRENT LIMITATION:    does not do hi-res design building; at TR
            onsets = round(onsets);

            if isempty(onsets), disp('No onsets found. Returning basis set only.'); end
            
            switch basistype
                case 'spm+disp'
                    bf = spm_get_bf(struct('dt', tr, 'name', 'hrf (with time and dispersion derivatives)', 'length', 20));
                    bf = bf.bf;
                case 'spmgamma2'
                    bf = spm_get_bf(struct('dt', tr, 'name', 'Gamma functions', 'length', 20, 'order', 2));
                    bf = bf.bf;
                case 'hrf'
                    bf = spm_get_bf(struct('dt', tr, 'name', 'hrf', 'length', 20));
                    bf = bf.bf;
                case 'gamma3'
                    bf = spm_get_bf(struct('dt', tr, 'name', 'hrf', 'length', 20));
                    bf = [bf.bf [0;bf.bf(1:end-1)] [bf.bf(2:end); 0]];
                case 'gamma4'
                    bf = spm_get_bf(struct('dt',tr,'name','hrf','length',20));
                    bf = [bf.bf [0;bf.bf(1:end-1)] [bf.bf(2:end); 0] [0;0;bf.bf(1:end-2)]];
                    
                case 'pain10'
                    hrf = spm_hrf(tr, [6 16 1 1 20 0 32]);
                    bf = conv([1 1 1 1 1], hrf);
                    bf = [bf [0; 0;bf(1:end-2)] ];

                    hrf = spm_hrf(tr, [6 16 .5 1 20 0 36]);
                    last = round(36 ./ tr);
                    bf = [hrf(1:last) bf(1:last,:)];
                    
                case 'painX'
                    hrf = spm_hrf(tr, [6 16 1 1 20 0 32]);
                    bf = conv(ones(1,round(stimlength/tr)), hrf);
                    bf = [bf [0; 0;bf(1:end-2)] ];

                    hrf = spm_hrf(tr, [6 16 .5 1 20 0 36]);
                    last = round(36 ./ tr);
                    bf = [hrf(1:last) bf(1:last,:)];
                    % % figure; b = randn(3,1000); f = bf * b; plot(f)
                    
                case 'flat10'
                    % FIR (flat basis fcns) for 10-s response in 3
                    % segments; only works for TR = 2!!
                    bf = [1 1 0 0 0; 0 0 1 0 0; 0 0 0 1 1]';
                    
                otherwise
                    error('unknown basis set.')
            end

            epochdur = 1;       % specify in s(TRs?), 0 is event
            if epochdur > 0
                % in 1/2 tr bins right now.
                for i =1:size(bf, 2)
                    bf2(:, i) = conv(bf(:, i), ones(round(epochdur./.5), 1));
                end
                bf = bf2;
                %figure;plot(bf)
            end

            [rowsbf, nbf] = size(bf);

            % define number of rows in X, if not entered.
            if isempty(rowsz), rowsz = max(onsets) + rowsbf; end

            z=zeros(rowsz, nbf);
            %one = ones(1, nbf);

            toolong = find(onsets > rowsz);
            if ~isempty(toolong), disp('Warning: Onsets after end of session.'); onsets(toolong) = []; end

            % make design matrix
            X = [];

            for i=1:length(onsets)
                xx = z(1:onsets(i)-1, :); % zeros before trial
                xx = [xx; bf];       % add basis functions
                xx =[xx; z(onsets(i)+rowsbf:end, :)];

                d = rowsz - size(xx, 1) ;
                if d < 0
                    xx = xx(1:rowsz, :);
                elseif d > 0
                    xx = padarray(xx, d, 0, 'post');
                end

                X = [X xx];
            end

            X(:, end+1) = 1;

            if isempty(onsets), X = [];  end
            
            varargout{1} = X;
            if nargout > 1, varargout{2} = pinv(X); end
            if nargout > 2, varargout{3} = bf; end

            if doplot
                create_figure('Trial beta plots', 1, 2); 
                imagesc(X); title('Design matrix'); colormap gray; drawnow;
            end
        otherwise
            error('Unknown output type.  Allowed: ''betas'', ''design''');
    end
end


function [h, t, w, auc] = get_beta_htw(b, doplot, bf)
    % need basis set (redundant code)
    %bf = spm_get_bf(struct('dt', .5, 'name', 'hrf (with time and
    %dispersion derivatives)', 'length', 20));
    %bf = bf.bf;

    [rowsbf, nbf] = size(bf);

    indx = 1;
    for i = 1:nbf:(length(b)-nbf+1)
        whbetas = i:i+nbf-1;
        ft = bf * b(whbetas);                    % fit for this trial
        [h(indx, 1), t(indx, 1), w(indx, 1), wt, hh, auc(indx, 1)] = fir2htw2(ft, length(ft), 0);     % height, time, and width for this trial
        indx = indx + 1;

    end

    if doplot
        create_figure('Trial beta plots', 1, 2, 1);
        subplot(1, 2, 2);
        
        plot(ft, 'k', 'LineWidth', 2);
        fir2htw2(ft, length(ft), 1);  % plot it
        title('Last Trial');
        drawnow
    end
end

