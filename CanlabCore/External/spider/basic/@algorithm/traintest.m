function [tstRes, trained, trnRes] = traintest(untrained, trnData, tstData, in4, in5)
%         [tstRes, algo, trnRes] = traintest(algo, trnData, tstData [, loss])
%   OR:   [tstRes, algo, trnRes] = traintest(algo, dat, trnInd, tstInd [, loss])
% 
% Trains algorithm algo on trnData [or on get(dat, trnInd)], and then tests the
% trained algorithm on tstData [or on get(dat, tstInd)].
% 
% For convenience, you may pass groups of data objects, or cell arrays of data objects,
% in place of single data objects. Groups of results and groups of trained algorithms
% are then returned as appropriate. (Note that this is less memory-efficient than looping
% through the train and/or test sets yourself.)
% 
% Other than that, the results are as you would expect from a call to algorithm/train
% followed by a call to algorithm/test.

if nargin < 5, in5 = []; end
if nargin < 4, in4 = []; end
if nargin < 3, tstData = []; end

if ~any(strcmp(class(trnData), {'data', 'data_global', 'group', 'cell'}))
	error('input argument trnData must be a data set, or a cell array/group of data sets')
end
if isempty(tstData), error('no test set specified'), end
if isa(tstData, 'double')
	index_mode = 1;
	trnIndex = tstData;
	if isempty(in4) | ~isa(in4, 'double'), error('no test indices supplied'), end
	tstIndex = in4;
	losstype = in5;
	tstData = [];
elseif any(strcmp(class(tstData), {'data', 'data_global', 'group', 'cell'}))
	index_mode = 0;
	losstype = in4;
	if ~isempty(in5), error('too many input arguments'), end
else
	error('input argument tstData must be a data set, or a cell array/group of data sets')
end

%%%% SGE support
if isdeferred(untrained)
	def = untrained.deferred;
	untrained.deferred = [];
	trained = untrained;
	def = reset(def, 'description', [mfilename ' ' get_name(untrained)]);
	if index_mode
		trnRes = get(trnData, trnIndex);
		tstRes = get(trnData, tstIndex);
		[tstRes.deferred trained.deferred trnRes.deferred] = qsub(def, max(nargout,1), mfilename, untrained, trnRes, tstRes, losstype);
	else
		trnRes = trnData;
		tstRes = tstData;
		[tstRes.deferred trained.deferred trnRes.deferred] = qsub(def, max(nargout,1), mfilename, untrained, trnData, trnIndex, tstIndex, losstype);
	end
	return
end
%%%%


% deal with the case that trnData and/or tstData are cell arrays of data objects, or groups of data objects
% (more convenient for the user in some cases, but slower and more memory-intensive)
if isa(trnData, 'group'), trnData = struct(trnData); trnData = trnData.child; end
if isa(tstData, 'group'), tstData = struct(tstData); tstData = tstData.child; end
if isa(trnData, 'cell') | isa(tstData, 'cell')
	if ~isa(trnData, 'cell'), trnData = {trnData}; end
	if ~isa(tstData, 'cell'), tstData = {tstData}; end
end
if isa(trnData, 'cell')
	if index_mode
		for i = 1:prod(size(trnData))
			[trnRes{i} trained{i}] = train(untrained, get(trnData{i}, trnIndex), losstype);
			tstRes{i} = test(trained{i}, get(trnData{i}, tstIndex), losstype);
		end
	else
		ntrn = prod(size(trnData));
		ntst = prod(size(tstData));
		if ntrn ~= 1 & ntst ~= 1 & ntrn ~= ntst, error('mismatched number of train & test datasets'), end
		for i = 1:ntrn
			[trnRes{i} trained{i}] = train(untrained, trnData{i}, losstype);
		end
		for i = 1:max(ntrn, ntst)
			itrn = min(i, ntrn);
			itst = min(i, ntst);
			tstRes{i} = test(trained{itrn}, tstData{itst}, losstype);
		end
	end
	if prod(size(trnRes)) == 1, trnRes = trnRes{1}; else trnRes = group(trnRes); end
	if prod(size(tstRes)) == 1, tstRes = tstRes{1}; else tstRes = group(tstRes); end
	if prod(size(trained)) == 1, trained = trained{1}; else trained = group(trained); end

else % deal with the single-data-set case

	if index_mode
		[trnRes trained] = train(untrained, get(trnData, trnIndex), losstype);
		tstRes = test(trained, get(trnData, tstIndex), losstype);
	else
		[trnRes trained] = train(untrained, trnData, losstype);
		tstRes = test(trained, tstData, losstype);
	end

end
