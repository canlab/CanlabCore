function out = extract(g, subfield, dim)
% EXTRACT(A, SUBFIELD [,DIM]) extract and group subfields from child objects
% 
% P = EXTRACT(ALGO, SUBFIELD)
% 
% If ALGO is an object that contains multiple 'child' objects (examples:
% GROUP, CHAIN, PARAM, CV), then this method is a way of extracting and
% grouping together some sub-property of all the child objects in one go.
% If the sub-property value is an ALGORITHM object for all the children,
% then the results are put into a GROUP. Otherwise, they are returned as a
% cell array.
% 
% Example:
%    C = cv(rfe(svm, 'feat=10; speed=15'));
%    [D C] = train(C, gen(toy));
%    SS = extract(C, 'child');
%    SVS = extract(C, 'child.Xsv');
%    R = extract(C, 'rank');
% 
%    % SS is now a GROUP of trained SVM objects, one per fold, that were
%    % wrapped in the RFE objects.
%    %
%    % SVS is now a GROUP of DATA objects, one per fold, containing the
%    % support vectors from the trained SVMs inside the RFE objects.
%    %
%    % R is now a cell array containing a feature ranking for each fold.
% 
% P = EXTRACT(ALGO, SUBFIELD, DIM)
% 
% This is the same, except that the resulting cells are concatenated in
% dimension DIM. In the above example,
%     R = extract(C, 'rank', 1);
%     ERR = extract(loss(D), 'Y', 2);
% would arrange the five rankings into the rows of a matrix R instead of
% keeping them in cells, and return the error rates in row vector ERR.

out = {};
all_algs = 1;
chn = subsref(g, struct('type', '.', 'subs', 'child'));
for i = 1:prod(size(chn))
	ch = chn{i};
	le = lasterr; lasterr('')
	s = subfield;
	err = '';
	while ~isempty(s)
		ind = min(find([s '.']=='.'));
		ss = s(1:ind-1);
		s = s(ind+1:end);
		eval('ch = subsref(ch, struct(''type'', ''.'', ''subs'', ss));', 'err=lasterr;')
		if ~isempty(err), error(['could not get subfield ''' subfield '''']), end
	end
	lasterr(le)
	out{end+1} = ch;
	all_algs = all_algs & isa(ch, 'algorithm');
end
if all_algs
	out = group(out);  % return as group
elseif nargin >= 3
	out = cat(dim, out{:}); % concatenate (numerics, structs, strings, etc)
else
	% leave as cell array
end
