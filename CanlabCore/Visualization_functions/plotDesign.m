function [X,d,out,handles] = plotDesign(ons,rt,TR,varargin)
% [X,d,out,handles] = plotDesign(ons,rt,TR,varargin)
%
% simple function to plot a design
% plots regressors and color-coded onset times as little sticks, with RT represented as height of the stick
%
% ons is a cell array of onset times in s   OR a delta indicator function matrix
% rt is a cell array of rts for each onset event
% TR is the repetition time for sampling, in s
% optional argument is the y offset for plotting rts, default = 2
%
% returns the model matrix (X) and the delta function d
%
% optional arguments
% 1     yoffset: default is 2
% 2     vector of epoch durations in sec for each trial type, default is events
%
% examples:
% plot epochs of different lengths stored in conditions(*).stimlength
% [X3,d] = plotDesign(evtonsets,[],1,2,cat(2,conditions.stimlength));

rtin = rt;      % original rt, to tell if rt is values or empty
yoffset = 2;
durs = []; out = [];
colors = {'r' 'g' 'b' 'c' 'm' 'y'};
samefig = 0;
basisset = 'hrf';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            % functional commands
            case 'yoffset', yoffset = varargin{i+1};
            case 'durs', durs = varargin{i+1}; durs = durs ./ TR;

            case {'color', 'colors'}, colors = varargin{i+1};
            case 'samefig', samefig = 1;
            case 'basisset', basisset =  varargin{i+1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% build models

if iscell(ons)
    if ~isempty(durs)
        for i = 1:length(ons)
            if length(durs) == 1
                xons{i} = [ons{i} repmat(durs, size(ons{i}, 1), 1)];
            end
        end
    else
        xons = ons; % xons is what we need to pass in to handle dur option
    end

    [X,d] = onsets2fmridesign(xons,TR, [], basisset);
else
    % note: will not convolve properly with durs
    if ~isempty(durs)
        warning('Must enter cell array of onsets if using dirs'); 
    end
    X = getPredictors(ons,spm_hrf(TR)./max(spm_hrf(TR)));
    d = ons; clear ons
    for i = 1:size(d,2)
        ons{i} = (find(d(:,i)) - 1) .* TR;
    end
end

if ~isempty(rtin),
    [X2,d2,out] = rt2delta(ons,rt,TR);
else
    % placeholder for plotting only
    for i = 1:length(ons), rt{i} = 1000 * ones(size(ons{i}));  end
end

while size(X,2) > length(colors), colors = [colors colors]; end


% make figure

if ~samefig, figure('Color','w'); end

if ~isempty(rtin), subplot(4,1,1); end
set(gca,'FontSize',16); hold on;
handles = [];

for i = 1:length(ons)

    if ischar(colors{i})
        handles(i) = plot(X(:,i),colors{i});
        h = plot([ons{i} ons{i}]'./(TR),[repmat(-yoffset,length(rt{i}),1) rt{i}]'./(1000.*TR) - yoffset,colors{i});

    else
        handles(i) = plot(X(:,i),'Color', colors{i});
        h = plot([ons{i} ons{i}]'./(TR),[repmat(-yoffset,length(rt{i}),1) rt{i}]'./(1000.*TR) - yoffset,'Color', colors{i});

    end


    if ~isempty(durs)
        for j = 1:length(ons{i})
            hh = drawbox(ons{i}(j)./TR,durs(i),colors{i},yoffset);
        end
    end

end

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
    if TR == 1, xlabel('Time (s)'), else xlabel('Time (TRs)'),end
end


return



function h1 = drawbox(time,dur,color,yoffset)

x = [0 1 1 0]; x = x * dur + time;
y = [0 0 1 1] - yoffset;

h1 = fill(x,y,color,'FaceAlpha',.5,'EdgeColor','none');

return

