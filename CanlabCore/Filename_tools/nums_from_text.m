function [nums,whnums] = nums_from_text(mytext)
% returns numeric values of numbers embedded in text strings
% given text string input
%
% :Usage:
% ::
%
%     function [nums,whnums] = nums_from_text(mytext)
%
% ..
%    tor wager
% ..


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

