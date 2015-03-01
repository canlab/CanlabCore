function x = equal(label1, label2)
% EQUAL(LABEL1, LABEL2) tests LABEL1 and LABEL2 for equality
% we need a special function since the two labels may be of different
% length and '==' might exit telling that the matrix dimensions don't agree
%
% returns 1 if true and 0 if false

% Authors: Alex J. Smola
% Created: 04/18/97
% Updated: 05/08/00
%
% This code is released under the GNU Public License
% Copyright by GMD FIRST and The Australian National University

if (length(label1) == length(label2))
	tmp = label1 == label2;
	x = (sum(tmp) == length(tmp));
else
	x = 0;
end


