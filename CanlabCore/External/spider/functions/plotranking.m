function plotranking(scores, X, dims)
% PLOTRANKING(DAT)
% PLOTRANKING(SCORES, X)
% 
% DAT is a data object returned by the MRANK method TEST.
% 
% SCORES may be one column [scores], two columns [indices, scores] or three
% columns [indices, scores, ranks] (the latter is the format returned by
% READRANKING or in the result field of an MRANK algorithm).
% 
% SCORES may also be a string, in which case it is automatically passed
% through READRANKING (to plot the contents of an output file).
% 
% Use NaN in the scores column to mark labelled points.
%
% X may be a string: if so, a variable X is loaded from the named file.
% If X is omitted, 'twomoons' is assumed. The first two columns of X
% provide the coordinates for the plotted data points.

if nargin < 3, dims = [1 2]; end
if nargin < 2, X = ''; end
if nargin < 1, scores = ''; end

if nargin == 1 & (isa(scores, 'data') | isa(scores, 'data_global'))
	X = get_x(scores);
	scores = get_y(scores);
end

if isempty(X), X = 'twomoons'; end
if isstr(X), load(X, 'X'); end
if isa(X, 'data') | isa(X, 'data_global'), X = get_x(X); end
X = X(:, dims);

if isstr(scores), scores = readranking(scores); end

m = size(X, 1);

if isempty(scores), scores = [1:m]'; end

if size(scores, 2) == 1
	if size(scores, 1) ~= m, error('SCORES must have one element per data point unless indices are supplied in the first column'), end
	scores = [[1:m]' scores];
end
if size(scores, 2) == 2
	labelled = isnan(scores(:, 2));
	scores(labelled, 3) = nan;
	t = compute_ranks(-scores(~labelled, 2));
	scores(~labelled, 3) = t(:);
end

X = X(scores(:, 1), :);
ranks = scores(:, 3);
scores = scores(:, 2);

labelled = isnan(scores);

plotpoints(X(~labelled, :), ranks(~labelled))
plot(X(labelled, 1), X(labelled, 2), 'marker', 'o', 'color', [1 0 0], 'markerfacecolor', [0 1 0], 'markersize', 10, 'linestyle', 'none')
axis square
figure(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotpoints(X, ranks)

fullcolour = 0;
m = size(X, 1);
cmap = jet(m);
[ranks i] = sort(max(ranks)-ranks);
X = X(i, :);
s = 10 + 90*ranks/m;
c = [1:m]';
if ~fullcolour, c = round(m/4)+round(m/2)*(c<m/2); end
c = cmap(c, :);
washeld = ishold;
plot(X(:, 1), X(:, 2), 'color', [1, 0,0])
hold on
if fullcolour
	scatter(X(:, 1), X(:, 2), 20, c, 'filled')
else
	scatter(X(:, 1), X(:, 2), s, c)
end
if ~washeld, hold off, end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rm = compute_ranks(am)

rm = zeros(size(am));
for j = 1:size(am, 2)
    a = am(:, j);
    [sorted r] = sort(a);
    r(r) = 1:length(r);
    while 1
        tied = sorted(min(find(diff(sorted) == 0)));
        if isempty(tied), break, end
        sorted(sorted==tied) = nan;
        r(a==tied) = mean(r(a==tied));
    end
    rm(:, j) = r;
end
