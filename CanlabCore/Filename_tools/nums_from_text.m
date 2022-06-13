function [nums,whnums] = nums_from_text(mytext)
% returns numeric values of numbers embedded in text strings
% given text string input
%
% :Usage:
% ::
%
%     function [nums,whnums] = nums_from_text(mytext)
%
%     see also regexp, which is much faster and very flexible.
%     e.g., to find numbers in text:
%     % Test string:
%     str = {'aix 1 2005' '99 db' '47_db' '1982 db 37'}'
%
%     % Convert matches to numbers:
%     S = regexp(str, '\d+', 'match', 'forceCellOutput');  S = [S{:}]'; nums = cellfun(@str2num, S);
%
%     % Take only the first match
%     S = regexp(str, '\d+', 'match', 'once', 'forceCellOutput');  nums = cellfun(@str2num, S);
%
%     % Say you have 2-digit or 4-digit years, and want to include only
%     % these in your match:
%     S = regexp(str, '\d{2,4}', 'match', 'once', 'forceCellOutput'); nums = cellfun(@str2num, S);
%
% ..
%    tor wager
% ..







mytext = strrep(mytext, ' ', '_');  % fix; did not work with spaces

nums = NaN;
ind = 1;
numind = 1;
mytext = [mytext ' '];

whnums = false(size(mytext));

for i = 1:length(mytext)

	tmp = str2num(mytext(i));

	if ~isempty(tmp) && isreal(tmp)
		mynum(ind) = tmp;
		catflag = 1;
		ind = ind + 1;

        whnums(i) = 1;
	else
		catflag = 0;
		ind = 1;

		if exist('mynum') == 1
			nums(numind) = catnums(mynum);
			numind = numind + 1;
		end

		clear mynum

	end
end

whomit = zeros(1,length(nums)); 
for k = 1:length(nums), if ~isreal(nums(k)), whomit(k) = 1; end, end
nums(find(whomit)) = [];
if isempty(nums), nums = NaN; end

return



function mynum = catnums(mynum)

	if length(mynum) > 1
		mult = 10 * ones(1,length(mynum));
		mult = mult .^ (length(mynum)-1:-1:0);
	else
		mult = 1;
	end

	mynum = sum(mynum .* mult);

return

