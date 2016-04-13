function hh = errorbar_horizontal(varargin)
% Error bar plot
%
% :Usage:
% ::
%
%    errorbar_horizontal(X,Y,L,U)
%
% plots the graph of vector X vs. vector Y with
% error bars specified by the vectors L and U.  L and U contain the
% lower and upper error ranges for each point in Y.  Each error bar
% is L(i) + U(i) long and is drawn a distance of U(i) above and L(i)
% below the points in (X,Y).  The vectors X,Y,L and U must all be
% the same length.  If X,Y,L and U are matrices then each column
% produces a separate line.
% ::
%
%    errorbar_horizontal(X,Y,E)
%    % or
%    errorbar_horizontal(Y,E)
%
% plots Y with error bars [Y-E Y+E].
% ::
%
%    ERRORBAR(...,'LineSpec')
%
% uses the color and linestyle specified by
% the string 'LineSpec'.  The color is applied to the data line and
% error bars while the linestyle and marker are applied to the data 
% line only.  See PLOT for possibilities.
% ::
%
%    errorbar_horizontal(AX,...)
%
% plots into AX instead of GCA.
% ::
%
%    H = errorbar_horizontal(...)
%
% returns a vector of errorbarseries handles in H.
%
% :Examples: To draws symmetric error bars of unit standard deviation
% ::
%
%    x = 1:10;
%    y = sin(x);
%    e = std(y)*ones(size(x));
%    errorbar(x,y,e)
%
% If means(:, 1) is x-values, means(:, 2) is y-values, stds(:, 1) is error
% on x, and stds(:, 2) is error on y, then:
% ::
%
%    lineh = errorbar_horizontal(means(:, 1), means(:, 2), stds(:, 1));
%    linehx = errorbar(means(:, 1), means(:, 2), stds(:, 2));
%
% ..
%    NOTE: Tor Wager modified from Matlab's ERRORBAR function to plot
%    HORIZONTAL error bars.
% ..

[v6,args] = usev6plotapi(varargin{:});
if v6
    warning(['MATLAB:', mfilename, ':DeprecatedV6Argument'],...
        ['The ''v6'' argument to %s is deprecated,',...
        ' and will no longer be supported in a future release.'], upper(mfilename));
  h = Lerrorbarv6(args{:});
else
  [cax,args,nargs] = axescheck(args{:});
  error(nargchk(1,inf,nargs,'struct'));
  [pvpairs,args,nargs,msg] = parseargs(args);
  if ~isempty(msg), error(msg); end %#ok
  error(nargchk(2,4,nargs,'struct'));

  hasXData = nargs ~= 2;
  x = [];
  switch nargs
   case 2
    [y,u] = deal(args{1:nargs});
    u = abs(u);
    l = u;
   case 3
    [x,y,u] = deal(args{1:nargs});
    if min(size(x))==1, x = x(:); end
    u = abs(u);
    l = u;
   case 4
    [x,y,l,u] = deal(args{1:nargs});
    if min(size(x))==1, x = x(:); end
  end
  if min(size(u))==1, u = u(:); end
  if min(size(l))==1, l = l(:); end
  if min(size(y))==1, y = y(:); end
  n = size(y,2);
  
  % Make sure that x,y,l and u all are the same size:
  if isempty(x)
      x = ones(size(y));
  end
  if ~isequal(size(x),size(y),size(l),size(u))
      error('MATLAB:errorbar:InputSizeMisMatch',...
          'X, Y and error bars must all be the same length');
  end

  % handle vectorized data sources and display names
  extrapairs = cell(n,0);
  if ~isempty(pvpairs) && (n > 1)
    [extrapairs, pvpairs] = vectorizepvpairs(pvpairs,n,...
                                            {'XDataSource','YDataSource',...
                        'UDataSource','LDataSource',...
                        'DisplayName'});
  end

  if isempty(cax) || isa(handle(cax),'hg.axes')
      cax = newplot(cax);
      parax = cax;
  else
      parax = cax;
      cax = ancestor(cax,'Axes');
  end

  h = []; 
  autoColor = ~any(strcmpi('color',pvpairs(1:2:end)));
  autoStyle = ~any(strcmpi('linestyle',pvpairs(1:2:end)));
  xdata = {};
  for k=1:n
    % extract data from vectorizing over columns
    if hasXData
      xdata = {'XData', double(full(x(:,k)))};  % tor modified to avoid private fcn
    end
    [ls,c,m] = nextstyle(cax,autoColor,autoStyle,k==1);
    
    h = [h torDrawHorizErrBars(y, x, l, u, cax)];
    
% %     h = [h specgraph.errorbarseries('YData',double(full(y(:,k))),...
% %                                     'UData',double(full(u(:,k))),...
% %                                     'LData',double(full(l(:,k))),xdata{:},...
% %                                     'Color',c,'LineStyle',ls,'Marker',m,...
% %                                     pvpairs{:},extrapairs{k,:},'parent',parax)];

  end
%   if autoColor
%     set(h,'CodeGenColorMode','auto');
%   end
%   set(h,'RefreshMode','auto');
  plotdoneevent(cax,h);
  h = double(h);
end

if nargout>0, hh = h; end

function h = Lerrorbarv6(varargin)
% Parse possible Axes input
error(nargchk(2,6,nargin,'struct'));
[cax,args,nargs] = axescheck(varargin{:});

x = args{1};
y = args{2};
if nargs > 2, l = args{3}; end
if nargs > 3, u = args{4}; end
if nargs > 4, symbol = args{5}; end

if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
    if nargs > 2,
        if ~ischar(l)
            l = l(:);
        end
        if nargs > 3
            if ~ischar(u)
                u = u(:);
            end
        end
    end
else
  npt = size(x,1);
end

if nargs == 3
    if ~ischar(l)  
        u = l;
        symbol = '-';
    else
        symbol = l;
        l = y;
        u = y;
        y = x;
        n = size(y,2);
        x(:) = (1:npt)'*ones(1,n);
    end
end

if nargs == 4
    if ischar(u),
        symbol = u;
        u = l;
    else
        symbol = '-';
    end
end


if nargs == 2
    l = y;
    u = y;
    y = x;
    n = size(y,2);
    x(:) = (1:npt)'*ones(1,n);
    symbol = '-';
end

u = abs(u);
l = abs(l);
    
if ischar(x) || ischar(y) || ischar(u) || ischar(l)
    error(id('NumericInputs'),'Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) || ~isequal(size(x),size(l)) || ~isequal(size(x),size(u)),
  error(id('InputSizeMismatch'),'The sizes of X, Y, L and U must be the same.');
end

tee = (max(x(:))-min(x(:)))/100;  % make tee .02 x-distance for error bars
xl = x - tee;
xr = x + tee;
ytop = y + u;
ybot = y - l;
n = size(y,2);

% Plot graph and bars
cax = newplot(cax);
hold_state = ishold(cax);

% build up nan-separated vector for bars
xb = zeros(npt*9,n);
xb(1:9:end,:) = x;
xb(2:9:end,:) = x;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xl;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = ytop;
yb(5:9:end,:) = ytop;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ybot;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;

[ls,col,mark,msg] = colstyle(symbol); 
if ~isempty(msg), error(msg); end %#ok
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

% ERRORBAR calls the 'v6' version of PLOT, and temporarily modifies global
% state by turning the MATLAB:plot:DeprecatedV6Argument warning off and
% on again.
oldWarn = warning('query','MATLAB:plot:DeprecatedV6Argument');
warning('off','MATLAB:plot:DeprecatedV6Argument');
try
    h = plot('v6',xb,yb,esymbol,'parent',cax); hold(cax,'on')
    h = [h;plot('v6',x,y,symbol,'parent',cax)]; 
catch
    warning(oldWarn);
    rethrow(lasterror);
end
warning(oldWarn);

if ~hold_state, hold(cax,'off'); end

function [pvpairs,args,nargs,msg] = parseargs(args)
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);
% check for LINESPEC
if ~isempty(pvpairs)
  [l,c,m,tmsg]=colstyle(pvpairs{1},'plot');
  if isempty(tmsg)
    pvpairs = pvpairs(2:end);
    if ~isempty(l) 
      pvpairs = {'LineStyle',l,pvpairs{:}};
    end
    if ~isempty(c)
      pvpairs = {'Color',c,pvpairs{:}};
    end
    if ~isempty(m)
      pvpairs = {'Marker',m,pvpairs{:}};
    end
  end
end
msg = []; % checkpvpairs(pvpairs); tor disabled; private fcn
nargs = length(args);

function str=id(str)
str = ['MATLAB:errorbar:' str];

function [l,c,m] = nextstyle(ax,autoColor,autoStyle,firsttime)
%NEXTSTYLE Get next plot linespec
%   [L,C,M] = NEXTSTYLE(AX) gets the next line style, color
%   and marker for plotting from the ColorOrder and LineStyleOrder
%   of axes AX.
%
%   See also PLOT, HOLD

%   [L,C,M] = NEXTSTYLE(AX,COLOR,STYLE,FIRST) gets the next line
%   style and color and increments the color index if COLOR is true
%   and the line style index if STYLE is true. If FIRST is true
%   then start the cycling from the start of the order unless HOLD
%   ALL is active.

if nargin == 1
  autoColor = true;
  autoStyle = true;
  firsttime = false;
end

co = get(ax,'ColorOrder');
lo = get(ax,'LineStyleOrder');

ci = [1 1];
if (isappdata(ax,'PlotHoldStyle') && getappdata(ax,'PlotHoldStyle')) || ...
      ~firsttime
  if isappdata(ax,'PlotColorIndex')
    ci(1) = getappdata(ax,'PlotColorIndex');
  end
  if isappdata(ax,'PlotLineStyleIndex')
    ci(2) = getappdata(ax,'PlotLineStyleIndex');
  end
end

cm = size(co,1);
lm = size(lo,1);

if isa(lo,'cell')
  [l,c,m] = colstyle(lo{mod(ci(2)-1,lm)+1});
else
  [l,c,m] = colstyle(lo(mod(ci(2)-1,lm)+1));
end
c = co(mod(ci(1)-1,cm)+1,:);

if autoStyle && (ci(1) == cm)
  ci(2) = mod(ci(2),lm) + 1;
end
if autoColor
  ci(1) = mod(ci(1),cm) + 1;
end
setappdata(ax,'PlotColorIndex',ci(1));
setappdata(ax,'PlotLineStyleIndex',ci(2));

if isempty(l) && ~isempty(m)
  l = 'none';
end
if ~isempty(l) && isempty(m)
  m = 'none';
end


function h = torDrawHorizErrBars(y, x, l, u, cax)
    
npt = size(y,1);

tee = (max(y(:))-min(y(:)))/100;  % make tee .02 y-distance for error bars
yl = y - tee;
yr = y + tee;

xtop = x + u;
xbot = x - l;
n = size(x,2);


% build up nan-separated vector for bars
yb = zeros(npt*9,n);
yb(1:9:end,:) = y;
yb(2:9:end,:) = y;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = yl;
yb(5:9:end,:) = yr;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = yl;
yb(8:9:end,:) = yr;
yb(9:9:end,:) = NaN;

xb = zeros(npt*9,n);
xb(1:9:end,:) = xtop;
xb(2:9:end,:) = xbot;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xtop;
xb(5:9:end,:) = xtop;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xbot;
xb(8:9:end,:) = xbot;
xb(9:9:end,:) = NaN;

[ls,col,mark,msg] = colstyle('-');
esymbol = ['-' col]; % Make sure bars are solid

h = plot(xb,yb,esymbol,'parent',cax); hold(cax,'on')

