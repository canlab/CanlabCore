function ff = prec_rec_plot(precisions, recalls, labels, varargin)
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
if nargin<1
    ff=prec_rec_plot([.92 .5 .1], [.984 .8 .05], ...
        {'Basophils', 'T-cells', 'B-cells'}', 'sizes', [42 24 83],...
        'colors', [1 1 0;.5 .5 0;0 0 1], 'invert', false);
    msg(Html.WrapHr(['<h2>Incorrect arguments</h2>'...
        'This is an example of a prec_rec_plot...<br><br>'...
        '(See the top of the file prec_rec_plot.m for this...)']),...
        0,'south east++', 'Incorrect arguments!');
    return
end
if nargin < 3 || isempty(labels)
    labels = 1:length(precisions);
    warning('Labels were not given for the precision-recall scatter plot!');
end

p=parseArguments();
parse(p,varargin{:});
args=p.Results;
invert = args.invert;
if args.visible
    visibility = 'on';
else
    visibility = 'off';
end
if max(args.sizes)>75 %normalize sizes for use with markers
    mx=max(args.sizes);
    args.sizes=25+(args.sizes/mx*50);
    nSizes=length(args.sizes);
    if nSizes>1 % sort biggest first for legend sake
        [~,I]=sort(args.sizes, 'descend');
        labels=labels(I);
        precisions=precisions(I);
        recalls=recalls(I);
        if size(args.colors,1)==nSizes
            args.colors=args.colors(I,:);
        end
        args.sizes=args.sizes(I);
    end
end
if invert
    xName = 'False positive rate';
    yName = 'False negative rate';
    precisions = 1 - precisions;
    recalls = 1 - recalls;
else
    xName = 'Precision';
    yName = 'Recall';
end

ff=figure('Name', [yName ' vs. ' xName],'visible', visibility);
ax=axes('Parent', ff);

try
    gscatter(ax, precisions, recalls, labels,args.colors,[],args.sizes);
catch ex 
    %r2019a and earlier do not support ax argument
    gscatter(precisions, recalls, labels,args.colors,[],args.sizes);
end
xlim(ax, [-.1 1.1]);
ylim(ax, [-.1 1.1]);
xlabel(ax, xName);
ylabel(ax, yName);
tickLabels={'\bf\color{blue}Perfect', ...
    '\bf\color[rgb]{.21 .21 .7}25%', '50%', ...
    '\bf\color{magenta}75%', ...
    '\bf\color{red}Failed'};
ticks=[0 .25 .5 .75 1];
if ~invert
    tickLabels=flip(tickLabels);
end
set(ax, 'xtick', ticks, ...
    'xTickLabel', tickLabels, 'xTickLabelRotation', -25, ...
    'ytick', ticks, ...
    'yTickLabel', tickLabels, 'yTickLabelRotation', -25);

function p=parseArguments()
        p = inputParser;
        addParameter(p,'invert',false,@(x) islogical(x));
        addParameter(p,'visible',true,@(x) islogical(x));
        addParameter(p, 'sizes', 25, @(x) all(x>=0) );
        addParameter(p, 'colors', {}, @(x) isnumeric(x) && size(x, 2)==3);
end
end