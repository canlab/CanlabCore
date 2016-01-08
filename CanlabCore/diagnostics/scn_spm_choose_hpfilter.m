function scn_spm_choose_hpfilter(spm_results_dir, varargin)
% Plots and choice of optimal high-pass filter from an SPM first-level
% model directory (with statistics and contrasts estimated.)
%
% :Usage:
% ::
%
%     scn_spm_choose_hpfilter(spm_results_dir, ['events_only'])
%
% SPM5 compatible and SPM8.
%
% Called by: scn_spm_design_check.m
% For all regressors or events only: see scn_spm_choose_hpfilter.m
%
% ..
%    Tor Wager
%    August 2010
% ..

if nargin < 1, spm_results_dir = pwd; end

spmfilename = fullfile(spm_results_dir, 'SPM.mat');
if ~exist(spmfilename, 'file')
    error('SPM.mat does not exist in %s\n', spm_results_dir);
end
load(spmfilename);

% Gets events of interest: All regressors, or events only if 'events_only'
% is input as keyword

wh_cols = scn_spm_get_events_of_interest(SPM, varargin{:});

% Prepare design columns of interest
X_interest = SPM.xX.X(:, wh_cols);
I = SPM.xX.X(:, SPM.xX.iB);
X_interest = X_interest - I * pinv(I) * X_interest;

k = size(X_interest, 2);

%% Start by plotting the regressors of interest

create_figure('regressors', 2, 2);
delete(subplot(2, 2, 1:2));
axh = axes('Position', [.05 .55 .9 .4]);
set(axh, 'FontSize', 16);

nregs = length(wh_cols);

if k < 20
    plot_matrix_cols(X_interest, 'horiz');
else
    imagesc(X_interest); colormap gray; colorbar
end


%% Plot the regressors after HP filtering in red

spm_hplength = SPM.xX.K(1).HParam;

X_filtered = filter_X(X_interest, cat(1, SPM.xX.K(:).X0));

hold on;
if k < 20
    han = plot_matrix_cols(X_filtered, 'horiz');
    set(han, 'Color', 'r');
else
    imagesc(X_filtered); colormap gray; colorbar
end



title(['Task regressors: Black = before filtering, Red = after current ' num2str(spm_hplength) ' s HP filter']);
drawnow

%% Calculate the cumulative power at different frequencies
TR = SPM.xY.RT;

[cumpower, myfreq, design_var_loss, five_percent_HPlength] = ...
    cumulative_power(X_interest, TR);

currenthp = mean(cat(1, SPM.xX.K.HParam));
best_hplen = five_percent_HPlength;

% Plot it
subplot(2, 2, 3);

for i = 1:nregs
    plot(myfreq(:, i), cumpower(:, i), 'k');
end
xlabel('Frequency (Hz)') ; ylabel('Cumpow (regs)');
axis tight
title(sprintf('<5%% regressor var. lost at HPlen = %3.0f sec', five_percent_HPlength))


try
    plot_vertical_line(currenthp, 'r');
    plot_vertical_line(1./best_hplen, 'b');
catch
    disp('Bug! check me')
    keyboard
end

drawnow

%% If Contrasts are entered, then plot those too
X_contrasts = get_contrast_vectors(SPM);
% can handle t- or F-contrasts

subplot(2, 2, 4);
if isempty(X_contrasts)
    title('NO Contrasts Entered!');
    
    
else
    [cumpower, myfreq, design_var_loss, five_percent_HPlength_contrasts] = ...
        cumulative_power(X_contrasts, TR);
    
    best_hplen = five_percent_HPlength_contrasts;
    
    % Plot it
    for i = 1:size(X_contrasts, 2)
        plot(myfreq(:, i), cumpower(:, i), 'k');
    end
    xlabel('Frequency (Hz)') ; ylabel('Cumpow (contrasts)');
    axis tight
    title(sprintf('<5%% con-variance lost at HPlen = %3.0f sec', five_percent_HPlength_contrasts))
    
    try
        plot_vertical_line(currenthp, 'r');
        plot_vertical_line(1./best_hplen, 'b');
    catch
        disp('Bug! check me')
        keyboard
    end
    
end

drawnow

%% Now re-filter the regressors with the best length

[Rmtx, dummy, K] = use_spm_filter(TR, size(X_interest, 1), 'none', 'specify', best_hplen);

X_filtered = filter_X(X_interest, K);

axes(axh); hold on;

if k < 20
    han = plot_matrix_cols(X_filtered, 'horiz');
    set(han, 'Color', 'b');
else
    imagesc(X_filtered); colormap gray; colorbar
end


text(1, nregs + .7, sprintf('Blue = filtering with <=5%% variance loss at %3.0f sec', best_hplen), 'FontSize', 16);
drawnow


end


function X_filtered = filter_X(X_interest, K)

%px = pinv(K);           % pinv of the filtering matrix
Rmtx = eye(size(K, 1)) - K * pinv(K);  % Residual-inducing matrix, y' = Rmtx*y

for i = 1:size(X_interest, 2)
    
    y = X_interest(:, i);             % select regressor to filter
    %y = y - K * px * y;     % residuals after filtering
    y = Rmtx * y;
    
    X_filtered(:, i) = scale(y);
    
end

end



function [cumpower, myfreq, design_var_loss, five_percent_HPlength] = cumulative_power(X_interest, TR)

for i = 1:size(X_interest, 2)
    
    
    [myfft, freq] = fft_calc(X_interest(:, i), TR);
    
    cumpower(:, i) = cumsum(myfft);
    myfreq(:, i) = freq;
    
    wh = find( cumpower(:, i) < .05 );
    if isempty(wh)
        wh = 2; % first one after intercept
        design_var_loss(i) = cumpower(wh, i);
        fprintf('No cutoff with < .5%% power for reg/contrast %3.0f! Min loss is %3.2f%%\n', i, 100*design_var_loss(i))
    end
    
    wh = wh(end);  % highest index value with < 5% power loss
    
    design_var_loss(i) = cumpower(wh, i);
    five_percent_loss_freq(i) = myfreq(wh, i);
    
end

% HP filter length that preserves 5% design variance loss
% or less in all regressors
five_percent_HPlength = 1 ./ min(five_percent_loss_freq);

end


function X_contrasts = get_contrast_vectors(SPM)
% can handle t- or F-contrasts
X_contrasts = [];
if ~isfield(SPM, 'xCon') || length(SPM.xCon) == 0
    return
end

% Remove intercepts
X = SPM.xX.X;
I = SPM.xX.X(:, SPM.xX.iB);
X = X - I * pinv(I) * X;

for mycon = 1:length(SPM.xCon)
    
    d = X * SPM.xCon(mycon).c;
    X_contrasts = [X_contrasts d];
    
end

end


