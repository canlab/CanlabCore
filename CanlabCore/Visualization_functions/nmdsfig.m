function f1 = nmdsfig(pc,varargin)
% ::
%
%    f1 = nmdsfig(pc,[opt. inputs in any order])
%
% Create a 1-D or 2-D plot with stimulus coordinates
%
% reserved keywords, each followed by appropriate input:
%   - case 'classes', clus = varargin{i+1};
%   - case 'names', names = varargin{i+1};
%   - case 'sig', sigmat = varargin{i+1};
%   - case 'thr', thr = varargin{i+1};
%   - case 'legend', legmat = varargin{i+1};
%   - case 'sig2', sigmat2 = varargin{i+1};a
%     sigmat2 can be thresholded at multiple values in thr
%   - case 'colors', colors = varargin{i+1};
%   - case 'sizes', sizes = varargin{i+1};
%   - case 'sigonly' plot regions with significant connections only
%   - case 'nolines', do not plot lines (lines plotted by default, but only if
%     sigmat is entered)
%
%     'linethickness', followed by matrix of line thickness values
%
%     NOTE: can enter sig matrix with non-zero values equal to line
%     thickness and use for both sig and linethickness inputs
%     % but thickness values should be scaled to integers for line thickness
%
% Creates a figure only if f1 output is requested
%
%  'fill', fill in areas around groups
%
%   **pc:**
%        is objects x dimensions
%
%   **clus:**
%        is a vector of object classes
%
%   **names:**
%        is cell array of names for rows of pc (objects), or empty ([])
%
%   **sig:**
%        is optional matrix of 1, -1, and 0 entries
%        signifies which pairs to connect with lines
%        positive elements are solid lines, negative elements are dashed
%          - Can be a series of t-maps in 3-D array
%
%   [opt] threshold vector of critical t-values, e.g., [2.2 5.4]
%   If used, enter t-maps in sigmat
%
%   [opt] a 2nd sigmat, if entered, will plot dashed lines instead of solid
%   ones. This is used by cluster_nmdsfig to plot interactions between
%   covariance and behavioral scores
%
% :Figure creation:
% If existing fig with tag 'nmdsfig', activates
% Otherwise, if fig handle requested as output, creates
% or if not, uses current figure.
%
% :Examples: c is output of cluster_nmdsfig
% ::
%
%    sizes = sum(c.STATS.sigmat);
%    f1 = nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',p_vs_c_heat.sig,'legend',{'Pos' 'Neg'},'sizes',sizes,'sizescale',[4 12]);
%
%    % add length legend
%    f1 =
%    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',p_vs_c_heat.sig,'legend',{'Pos' 'Neg'}, ...
%    'sizes',sizes,'sizescale',[4 12],'lengthlegend',c.r);
%
%    % Auto size scaling based on number of connections:
%    f1 =
%    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',c.STATS.sigmat,'legend',{'Pos' 'Neg'},'sizescale',[4 12],'lengthlegend',c.r);
%
%    f1 =
%    nmdsfig(c.GroupSpace,'classes',c.ClusterSolution.classes,'names',c.names,'sig',p_vs_c_heat.sig,'legend',{'Pos' 'Neg'},'sizescale',[4 16],'sigonly');
%
% :SEE ALSO: cluster_nmdsfig

% ..
%    set up defaults
% ..
pc = double(pc);

nobj = size(pc,1);
clus = ones(1,nobj);
for i = 1:size(pc,1), names{i} = ['V' num2str(i)]; end
sigmat = [];
sigmat2 = [];
thr = 0;
dolines = 1;


linecolor(1,:) = [0 0 0];   % first line color, 'positive'
linecolor(2,:) = [0 .3 1];   % second line color, 'negative'

linestyle{1} = '-';         % first line style, positive
linestyle{2} = '-';         % second line style, negative


legmat={'positive','negative'};
colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
sizes = 6 * ones(nobj,1);
dosizescale = 0;
dolengthlegend = 0;
selectpoints = 0;
dofill = 0;
dolinescale = 0;

for i = 1:length(varargin)
    if isstr(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'classes', clus = varargin{i+1};
            case 'names', names = varargin{i+1}; varargin{i + 1} = [];
            case 'sig', sigmat = varargin{i+1};
            case 'thr', thr = varargin{i+1};
            case 'legend', legmat = varargin{i+1};
            case 'sig2', sigmat2 = varargin{i+1};
            case 'colors', colors = varargin{i+1};
            case {'size', 'sizes'}, sizes = varargin{i+1};
            case 'sizescale', dosizescale = 1; minmax = varargin{i+1};
            case 'lengthlegend', dolengthlegend = 1; rmatx = varargin{i+1};
            case 'sigonly', selectpoints = 1;
            case 'nolines', dolines = 0;
            case 'fill', dofill = 1;

            case {'linescale', 'linethickness'}, dolinescale = 1; thicknessvals = varargin{i + 1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

thr = sort(thr);        % thresholds in ascending order, dark to light

% color  stuff
while max(clus) > length(colors), colors = [colors colors]; end
if length(colors{1}) == 2
    colormode = 'text';
elseif length(colors{1}) == 3
    colormode = 'vector';
else
    error('I don''t understand the format of the colors input.');
end

% size stuff: sigmoid scaling
%sizes = sum(c.STATS.sigmat);
if length(sizes) == 1
    sizes = repmat(sizes, nobj, 1);
end

if dosizescale
    % if we haven't entered sizes but we want to scale, assume we want to
    % use sig. connections
    if ~isempty(sigmat)
        sizes = sum(abs(sigmat));
    end
    if ~isempty(sigmat2)        % sig2 connections count for 1/2, remove already sig ones
        sizes = sizes + .5 * sum(abs(sigmat2 .* ~sigmat));
    end
    rg = abs(diff(minmax));
    zsz = zscore(sizes);
    sizes = min(minmax) + rg * (1 ./ (1+exp(-3*zsz)));
end

% --------------------------------------
%%% remove extra dimensions for plotting
% --------------------------------------
if size(pc,2) > 2
    disp('Using only first 2 dimensions of pc for plotting');
elseif size(pc,2)==1;
    disp('switching to 1D configuration plot');
    nmdsfig1D(pc,clus,names,varargin);
    return
end
pc=pc(:,1:2);


% --------------------------------------
%%% draw figure and set tag
% --------------------------------------
% Figure creation:
% If existing fig with tag 'nmdsfig', activates
% Otherwise, if fig handle requested as output, creates
% or if not, uses current figure.

f1 = findobj('Tag', 'nmdsfig'); % look for existing
if ~isempty(f1) && strcmp(get(f1,'Type'), 'figure')
    % we have one already, just clear
    figure(f1);
elseif nargout > 0
    f1 = create_figure;
else f1 = gcf;
end
set(f1,'Tag','nmdsfig');


% --------------------------------------
%%% fill areas, if requested
% --------------------------------------
if dofill
    areaHandles = nmdsfig_fill(clus, pc,  colors);
end



% --------------------------------------
% set up colour and line thickness parameters
% --------------------------------------


% % % if length(thr) > 1,
% % %     sigmat = repmat(sigmat,[1 1 length(thr)]);
% % %     % choose colors
% % %     sigcol = [.5 0];  %.5:-.5./length(thr):0;
% % % else
% % %     sigcol = 0;       % black, width = 3
% % % end

% --------------------------------------
% if multiple thr are entered, treats sigmat as a series of t-value maps
% replicate sigmat and threshold.
% draw AFTER 2ND THRESHOLD
% --------------------------------------
if ~isempty(sigmat) && dolines
    
    sigmat = double(sigmat);
    
    for i = 1:length(thr)
        wh = sign(sigmat(:,:,i)) .* (abs(sigmat(:,:,i)) >= thr(i));   % preserve pos or neg t-values
        
        if i == 1, sigfirstthr = wh; end

        line_handles = nmdsfig_tools('drawlines',pc, wh, linecolor, linestyle);
        
        hh1 = line_handles.hhp;
        hh2 = line_handles.hhn;

        %[hh1,hh2] = drawlines(pc,wh,sigcol(i),legmat);
    end
    
    
    if dolinescale
        pindx = 1;
        nindx = 1;
        
        for i = 1 : size(pc,1)
            for j = i+1 : size(pc,1)
                
                if wh(i,j) > 0
                    set(hh1(pindx), 'LineWidth', thicknessvals(i, j));
                    pindx = pindx + 1;
                    
                elseif wh(i,j) < 0
                    set(hh2(nindx), 'LineWidth', thicknessvals(i, j));
                    nindx = nindx + 1;
                    
                end
                
            end
        end
    end
            
else
    hh1 = []; hh2 = [];
end


% --------------------------------------
% if a 2nd sigmat is entered, plot dotted lines
%
% --------------------------------------
if ~isempty(sigmat2) && dolines
    
    sigmat2 = double(sigmat2);
    
    for i = 1:length(thr)
        wh = sign(sigmat2(:,:,i)) .* (abs(sigmat2(:,:,i)) >= thr(i));   % preserve pos or neg t-values
        wh = wh .* ~(logical(sigfirstthr));    % omit sig in first sigmat, first thresh

        line_handles = nmdsfig_tools('drawlines',pc, wh, linecolor, [{':'} {':'}]);
        hh3 = line_handles.hhp;
        hh4 = line_handles.hhn;

        %[hh3,hh4] = drawlines(pc,wh,sigcol(i),legmat,[0 0 0; 0 .5 1],[{':'} {':'}]);
    end
    
    legmat = [legmat legmat];
    hh = [hh1 hh2 hh3 hh4];
    
end

% --------------------------------------
% if multiple thr are entered, treats sigmat as a series of t-value maps
% replicate sigmat and threshold, then draw
% --------------------------------------
if ~isempty(sigmat) && dolines
    sigmat = double(sigmat);
    for i = 1:length(thr)
        wh = sign(sigmat(:,:,i)) .* (abs(sigmat(:,:,i)) >= thr(i));   % preserve pos or neg t-values

        line_handles = nmdsfig_tools('drawlines',pc,wh, linecolor, linestyle);
        hh1 = line_handles.hhp;
        hh2 = line_handles.hhn;
        %[hh1,hh2] = drawlines(pc,wh,sigcol(i),legmat);
    end
else
    hh1 = []; hh2 = [];
end

% do legend for 2nd
if ~isempty(sigmat2) && dolines
    wh = zeros(1,4);
    if isempty(hh1), wh(1) = 1; end
    if isempty(hh2), wh(2) = 1; end
    if isempty(hh3), wh(3) = 1; end
    if isempty(hh4), wh(4) = 1; end
    legmat(find(wh)) = [];
    legend(hh,legmat);
end

% --------------------------------------
% plot points
% --------------------------------------
plotthispoint = ones(nobj,1);
if selectpoints
    plotthispoint = (sum(abs(sigmat)) | sum(abs(sigmat2)));
end

hold on
for j = 1:nobj
    if plotthispoint(j)
        switch colormode
            case 'text'
                plot(pc(j,1),pc(j,2),colors{clus(j)},'MarkerFaceColor',colors{clus(j)}(1),'LineWidth',2,'MarkerSize' ,sizes(j))
            case 'vector'
                plot(pc(j,1),pc(j,2),'o','Color',colors{clus(j)},'MarkerFaceColor',colors{clus(j)},'LineWidth',2,'MarkerSize' ,sizes(j))
        end
    end
end
if ~isempty(names)
    for i =  1:nobj
        if plotthispoint(i)
            text(pc(i,1) + .025 * [max(pc(:,1))-min(pc(:,1))],pc(i,2)+.0,names{i},'Color','b','FontSize',18);
        end
    end
end

xlabel('Component 1'),ylabel('Component 2')
%set(gca,'XTickLabel',{[]});
%set(gca,'YTickLabel',{[]});

axis equal;
axis image;

if dolengthlegend
    nmdsfig_legend(pc,rmatx);
end

if dofill
    set(findobj(f1, 'Type', 'Text'),'Color','k','FontWeight','b')
end

return

% NOW IN nmdsfig_tools('drawlines')
%
% % % function [hhp,hhn] = drawlines(pc,sigmat,sigcol,legmat,varargin)
% % % % function [hhp,hhn] = drawlines(pc,sigmat,sigcol,legmat,[color],[style])
% % % %
% % % % pc is stim coords
% % % % sigmat is matrix of which lines to draw
% % % % sigcol is multiplier to scale line width (0 is fixed width)
% % % %   sigcol is a scalar between zero and one, lower is 'more salient'
% % % % legmat is not used now (for legend stuff)
% % % %
% % % % Example:
% % % % [hh3,hh4] = drawlines(pc,wh,sigcol(i),legmat,[.5 .5 .5; 0 1 1],[{':'}
% % % % {':'}]);
% % %
% % % % sigcol is a scalar between zero and one, lower is 'more salient'
% % % hold on
% % %
% % % % linew
% % % lw = 2;
% % %
% % % color(1,:) = [0 0 0];   % first line color, 'positive'
% % % color(2,:) = [0 .5 1];   % second line color, 'negative'
% % % style{1} = '-';         % first line style
% % % style{2} = '-';         % second line style
% % %
% % % if nargin < 3, sigcol = 0; end
% % % if length(varargin) > 0, color = varargin{1}; end
% % % if length(varargin) > 1, style = varargin{2}; sigcol = 0; end
% % %
% % % hhp=[];hhn=[];
% % % for i=1:size(pc,1);
% % %     for j=1:size(pc,1);
% % %         if sigmat(i,j) > 0
% % %             hhp(end+1) = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
% % %             set(hhp,'Color',color(1,:) + [1 1 1] * abs(sigcol),'LineStyle',style{1},'LineWidth',lw - (lw*sigcol)-1)
% % %         elseif sigmat(i,j) < 0
% % %             hhn(end+1) = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
% % %             set(hhn,'Color',color(2,:) + abs([sigcol sigcol 0]),'LineStyle',style{2},'LineWidth',lw - (lw*sigcol)-1)
% % %         end
% % %     end
% % % end
% % % %legend ([hhp hhn],legmat);
% % %
% % % % store in gui data for later deletion, etc.
% % % figh = gcf;
% % % linehandles = [hhp hhn];
% % % data = guidata(figh);
% % % if ~isfield(data,'linehandles')
% % %     data.linehandles = linehandles;
% % % else
% % %     data.linehandles = [data.linehandles linehandles];
% % % end
% % % guidata(figh,data);
% % %
% % % return



function f1 = create_figure
scnsize = get(0,'ScreenSize');
f1 = figure('position',round([50 50 scnsize(3)./2 scnsize(3)./2]),'color','white');
return
