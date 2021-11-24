function [X,d,out,handles] = plotDesign(ons,rt,TR,varargin)
% :Usage:
% ::
%
%    [X,d,out,handles] = plotDesign(ons,rt,TR,varargin)
%
% simple function to plot a design
% plots regressors and color-coded onset times as little sticks, with RT represented as height of the stick
%
% :Inputs:
%
%   **ons:**
%        a cell array of onset times in s
%
%        OR a delta indicator function matrix
%        Event durations (durs) will ONLY work with cell array inputs
%
%        * optional *
%        The second column of each cell of ons can be a series of event durations for
%        each event.
%
%   **rt:**
%        is a cell array of rts or other parametric modulator for each onset event (or empty if no values)
%
%   **TR:**
%        repetition time for sampling, in s
%
% :Optional Inputs:
% returns the model matrix (X) and the delta function d
%
%   optional arguments
%     1. y offset for plotting rts, default = 2
%     2. vector of epoch durations in sec for each trial type, default is events
%
%   **'yoffset':**
%        followed by yoffset for plotting rts; default is auto scale
%
%   **'durs':**
%        followed by durations in sec, either:
%
%        Constant duration
%
%        Vector of one duration for each event type
%        Cell array of one duration per trial
%        * Note: You can also add duration to ons input
%        instead; see ons above for more info *
%
%   **{'color', 'colors'}:**
%        followed by cell array of colors
%
%   **'samefig':**
%        keep on same figure
%
%   **'basisset':**
%        followed by name of basis set
%
%   **'overlapping':**
%        Default is to plot separate lines in separate
%        vertical positions. To plot overlapping in same location, enter this.
%
% :Examples: plot epochs of different lengths stored in conditions(*).stimlength
% ::
%
%    [X3,d] = plotDesign(evtonsets,[],1,2,cat(2,conditions.stimlength));
%
% ------------------------------------------------------------------------
% Create a single-trial design with 5 trials, TR = 1, 
% modeling each event with a 3-sec epoch
% ons = [{1} {5} {9} {13} {17} {21}];
% [X,d,out,handles] = plotDesign(ons,[], 1, 'durs', 3);
% 
% % Same as above, but with 2 sec epoch
% [X,d,out,handles] = plotDesign(ons,[], 1, 'durs', 2);
% 
% % Same as above, but with event-related (no epoch)
% [X,d,out,handles] = plotDesign(ons,[], 1);
% 
% % Same as above, but with TR = 1.3
% [X,d,out,handles] = plotDesign(ons,[], 1.3);
%
% ------------------------------------------------------------------------
% % Create random event-related design with an event every 3 sec, 4
% % conditions with 20% frequency each, 2 sec epochs, and plot it with 1.3 sec TR:
% ons = create_random_onsets(100, 3, [.2 .2 .2 .2], 2);
% [X,d,out,handles] = plotDesign(ons,[], 1.3);
%
% % Plot VIFs for this design:
% create_figure('vifs'); getvif(X, 0, 'plot');
%
% % Get efficiency
% nconditions = length(ons);
% 
% contrasts = create_orthogonal_contrast_set(nconditions);
% contrasts(:, end+1) = 0; % for intercept
% 
% e = calcEfficiency(ones(1, size(contrasts, 1)), contrasts, pinv(X), []);
%
% ------------------------------------------------------------------------
% % Create a single-trial design with an event every 10 sec
% % 4-sec epoch dur, and 1.3 sec TR.  plot it.
% ons = mat2cell([1:10:200]', ones(20, 1))';
% [X,d,out,handles] = plotDesign(ons,[], 1.3, 'durs', 4);
% 
% % Plot VIFs for this design:
% create_figure('vifs'); getvif(X, 0, 'plot');
% 
% % Get efficiency
% nconditions = length(ons);
% 
% contrasts = create_orthogonal_contrast_set(nconditions);
% contrasts(:, end+1) = 0; % for intercept
% 
% e = calcEfficiency(ones(1, size(contrasts, 1)), contrasts, pinv(X), []);
% ------------------------------------------------------------------------
%
% :See Also: onsets2fmridesign, create_random_er_design,
% create_block_design, create_random_onsets
%
% ..
%    Programmers' notes
%    DURS still needs some work to match onsets2fmridesign
%    tor edited 8/2015 to fix some ease-of-usage issues and document
% ..

rtin = rt;      % original rt, to tell if rt is values or empty
yoffset = [];
durs = [];
out = [];

% default colors
%colors = {'r' 'g' 'b' 'c' 'm' 'y'};
colors = get(gcf, 'DefaultAxesColorOrder');
colors = mat2cell(colors, ones(size(colors, 1), 1), 3);

samefig = 0;
basisset = 'hrf';

doseparatelines = 1;

inputs_to_pass = {}; % pass to onsets2fmridesign

% Process optional inputs
% ------------------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'yoffset', yoffset = varargin{i+1};
            case 'durs', durs = varargin{i+1};
                durs = durs ./ TR;
                if length(durs) == 1, durs = repmat(durs, 1, length(ons)); end
                
            case {'color', 'colors'}, colors = varargin{i+1};
            case 'samefig', samefig = 1;
            case 'basisset', basisset =  varargin{i+1};
                
            case 'nonlinsaturation', inputs_to_pass = {'nonlinsaturation'};
                
            case 'singletrial', inputs_to_pass = {'singletrial'};
                
            case 'overlapping', doseparatelines = 0;
                
            case 'nononlin' % ignore
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% Build models
% Create xons, onsets with durations added if needed
% ------------------------------------------------------------------------

if iscell(ons)
    
    [xons, ons_includes_durations] = parse_cell_onsets_durs(ons, durs);
    
    len = max(cat(1, xons{:}));
    
    len = len + 16; % add 16 sec past last event
    
    if length(len) > 1
        % in case durations are entered
        len = len(1) + len(2);
    end
    % length must be even multiple of TR
    if mod(len, TR)
        len = len + (TR - mod(len, TR));
        len = round(len * 100) ./ 100; % avoid precision errors later
    end
    
    [X,d] = onsets2fmridesign(xons, TR, len, basisset, inputs_to_pass{:});
    
else
    
    error('Enter onsets in sec in a cell array, one cell per event type.');
    
end

% RT model - special
% ------------------------------------------------------------------------

if ~isempty(rtin)
    [~, ~, out] = rt2delta(ons, rt, TR);
else
    % placeholder for plotting only
    % these are the amplitudes of the stick functions in the plot
    for i = 1:length(ons), rt{i} = 1000 * ones(size(ons{i}, 1), 1);  end
end


% make figure
% ------------------------------------------------------------------------

while size(X, 2) > length(colors), colors = [colors colors]; end

if ~samefig, figure('Color','w'); end

if ~isempty(rtin), subplot(4, 1, 1); end

set(gca,'FontSize',16); hold on;
handles = [];


% Plot lines for regressors
% ---------------------------------------------------------
% Plot in seconds (convert X from TRs)
time_in_sec = TR .* [0:size(X, 1)-1]'; %#ok<NBRAK>

if doseparatelines
    
    handles = plot_matrix_cols(X, 'horiz', time_in_sec, colors, 2);
    
else
    % overlapping on same plot
    for i = 1:length(ons)
        
        handles(i) = plot(time_in_sec, X(:,i), 'Color', colors{i});
        
    end
    
end


% Set limits for boxes and sticks
% ---------------------------------------------------------
%ymax = .1 * size(X, 2);  %.05 * (max(X(:)) - min(X(:)));
% ymax = .2; %* (max(X(:)) - min(X(:)));

if doseparatelines
    
    for i = 1:length(handles)
        yrange(i) = range(handles(i).YData);
    end
    
    ymax = .2 .* max(yrange);
    
else
    ymax = .2;
end

if isempty(yoffset)
    % default - auto
    yoffset = min(X(:)) - ymax;
end

for i = 1:length(xons)
    
    % for sticks
    xvals = [xons{i}(:, 1) xons{i}(:, 1)]'; %./ TR; % was ons
    
    if doseparatelines
        
        yvals1 = (i - 1) + ones(1, length(rt{i})) .* (max(handles(1).YData) + ymax);
        yvals2 = yvals1 - ymax .* (rt{i} ./ 1000)' ;
        yvals = [yvals1; yvals2];
        
    else
        
        yvals1 = ones(1, length(rt{i})) .* yoffset + ymax;
        yvals2 = yvals1 - ymax .* (rt{i} ./ 1000)' ;
        yvals = [yvals1; yvals2];
        
    end
    
    % onsets
    h = plot(xvals, yvals, 'Color', colors{i}, 'LineWidth', 3);
    
    
    % boxes for event durations
    
    if ons_includes_durations
        
        for j = 1:size(xons{i}, 1)
            
            yoffset = yvals2(j);
            
            % now durs are always integrated into ons above.
            hh = drawbox(xons{i}(j, 1), xons{i}(j, 2), yoffset, ymax, colors{i});
            
            set(hh, 'EdgeColor', 'none');
        end
    end % boxes
    
end % onsets

axis tight



title('Predicted activity')


% plot 2

if ~isempty(rtin)
    
    subplot(4,1,2); set(gca,'FontSize',16); hold on;
    for i = 1:length(ons)
        
        plot(out.rtlinearX(:,i),colors{i})
        h = plot([ons{i} ons{i}]'./(TR),[repmat(-yoffset,length(rt{i}),1) rt{i}]'./(1000.*TR) - yoffset,colors{i});
        
    end
    title('Activity x reaction time (linear)')
    
    subplot(4,1,3); set(gca,'FontSize',16); hold on;
    for i = 1:length(ons)
        
        plot(out.rtquadX(:,i),colors{i})
        h = plot([ons{i} ons{i}]'./(TR),[repmat(-yoffset,length(rt{i}),1) rt{i}]'./(1000.*TR) - yoffset,colors{i});
        
    end
    title('Activity x reaction time (quadratic)')
    
    subplot(4,1,4); set(gca,'FontSize',16); hold on;
    ind = 1;
    for i = 1:length(ons)
        
        hh(1) = plot(out.rtclassX(:,ind),colors{i},'LineStyle','-'); ind = ind + 1;
        hh(2) = plot(out.rtclassX(:,ind),colors{i},'LineStyle','--'); ind = ind + 1;
        hh(3) = plot(out.rtclassX(:,ind),colors{i},'LineStyle',':'); ind = ind + 1;
        
        %tmp = find(out.rtclass(:,ind));
        %h = plot([tmp tmp]'./(TR),[repmat(-yoffset,length(rt{i}),1) rt{i}]'./(1000.*TR) - yoffset,colors{i});
        
    end
    
    legend(hh,{'Fast' 'Medium' 'Slow'})
    xlabel('Time (TRs)')
    title('Activity for trials classified by RT')
    
else
    % single plot, add x label
    xlabel('Time (s)')
end


end % function



%
% function h1 = drawbox(time,dur,color,yoffset)
%
% x = [0 1 1 0]; x = x * dur + time;
% y = [0 0 1 1] - yoffset;
%
% h1 = fill(x,y,color,'FaceAlpha',.5,'EdgeColor','none');
%
% return



function [xons, ons_includes_durations] = parse_cell_onsets_durs(ons, durs)

xons = [];


[~, ons_includes_durations] = check_onsets(ons);

% sizes of onsets in each cell: First col is # events, 2nd is durations
% if they exist.
% sz = cellfun(@size, ons, 'UniformOutput', false)';
% sz = cat(1, sz{:});

% ons_includes_durations = any(sz(:, 2) > 1);  %  was: size(sz, 2) > 1;  sz(:, 2) = number of columns

if ons_includes_durations
    % we are done.
    xons = ons;
    
    if ~isempty(durs), warning('Onsets input already contains durations. Durs input will be ignored'); end
    
elseif ~isempty(durs) && iscell(durs)
    % We have durations for each event
    
    for i = 1:length(durs)
        xons{i} = [ons{i} durs{i}];
    end
    
    % update
    ons_includes_durations = true;
    
elseif ~isempty(durs) && length(durs) == 1
    % We have the same duration for every event
    
    for i = 1:length(ons)
        xons{i} = [ons{i} repmat(durs, size(ons{i}, 1), 1)];
    end
    
    % update
    ons_includes_durations = true;
    
elseif ~isempty(durs) && length(durs) == length(ons)
    % We have one duration value for each event type
    
    for i = 1:length(ons)
        xons{i} = [ons{i} repmat(durs(i), size(ons{i}, 1), 1)];
    end
    
    % update
    ons_includes_durations = true;
    
elseif ~isempty(durs)
    warning('Durs input is entered, but format/length is unrecognized. Will not be used.');
    
elseif isempty(durs)
    % No durations
    xons = ons; % xons is what we need to pass in to onsets2fmridesign handle dur option
    
else
    error('impossible case? Check code.');
    
end  % Parse what to do with durations

end % function


function [sz, ons_includes_durations] = check_onsets(ons)
% [sz, ons_includes_durations] = check_onsets(ons)

% sizes of onsets in each cell: First col is # events, 2nd is durations
% if they exist.
sz = cellfun(@size, ons, 'UniformOutput', false)';
sz = cat(1, sz{:});

ons_includes_durations = any(sz(:, 2) > 1);

end
