function scores = readranking(str)
% SCORES = READRANKING(STR)
% 
% Read ranking results from the string STR (or contents of file named STR),
% converting from the output format of the ranking program to a
% three-column matrix: indices, scores, ranks. Labelled points have score
% and rank set to NaN.
% 
% SCORES = READRANKING(SCORES)
% 
% Can also be used to convert a vector of scores into the three-column
% format.

s = [];
if nargin < 1, str = ''; end

if isnumeric(str)
	scores = str;
else
	if isempty(str)
		[a b] = uigetfile('*', 'Text file to open:');
		if ~isempty(a) & isstr(a) & isstr(b)
			str = fullfile(b, a);
		end
	end
	if isempty(str), return, end
	
	if ~any(str == 13 | str == 10)
		fid = fopen(str, 'rt');
		if fid == -1, error(['could not open ''' str '''']), end
		s = fscanf(fid, '%c', inf);
		s = char(s(:)');
		if fclose(fid) ~= 0, warning(['unable to close file ''' str '''']), end
		str = s;
	end
	str = sprintf('scores = [\n%s\n];', strrep(str, ':', ','));
	str = strrep(str, '#', sprintf('\n'));
	eval(str)
end

if size(scores, 2) == 1, scores = [[1:length(scores)]' scores]; end
labelled = isnan(scores(:, 2));
scores(labelled, 3) = nan;
t = compute_ranks(-scores(~labelled, 2));
scores(~labelled, 3) = t(:);


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
