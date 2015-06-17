function outlist = insert_rests(inlist,restevery,restlength,varargin)
% function outlist = insert_rests(inlist,restevery,restlength,numstim [if using variable restevery.])
%
% if restevery is variable, e.g. [7 8 9; 7 7 8 ;...], then insert them at fixed intervals.
%	for variable, rests in a single list are in a column.  columns index lists.

if size(restevery,1) > 1, type = 'variable';, else type = 'constant';,end

inlist = double(inlist);
restevery = double(restevery);

switch type
case 'variable'
% ========================================================================================================
	numstim = varargin{1};
	for i = 1:size(inlist,2)	% for each list
		index = 1; outindex = 1;
		list = inlist(:,i);
		for j = 1:size(restevery,1)
			if index+restevery(j,i)-1 > size(list,1), 
				tempoutlist(outindex:outindex+restevery(j,i)-1,1) = zeros(restevery(j,i),1);			% insert zeros if rest list is too long.
			else
				tempoutlist(outindex:outindex+restevery(j,i)-1,1) = list(index:index+restevery(j,i)-1);		% copy stimuli over
			end
			index = index + restevery(j,i); outindex = outindex + restevery(j,i);						% update indices
			tempoutlist(outindex:outindex+restlength-1,1) = 0;  											% insert rests
			outindex = outindex + restlength;															% update outlist index after rests
		end
		
		% adjust length, if necessary
		if size(tempoutlist,1) < numstim, tempoutlist = [tempoutlist; inlist(index:index+numstim-size(tempoutlist,1),1)];,end
		if size(tempoutlist,1) > numstim, tempoutlist = tempoutlist(1:numstim,1);,end						
		
		outlist(:,i) = tempoutlist;
	end

	
	
case 'constant'
% ========================================================================================================
outindex = restevery + 1;
for i = restevery+1:restevery:size(inlist,1)
    outlist(outindex-restevery:outindex-1,:) = inlist(i-restevery:i-1,:);
    outlist(outindex:outindex+restlength-1,:) = 0;  % insert rests
    outindex = outindex + restevery + restlength;
end

outlist(size(outlist,1)+1:size(outlist,1)+(size(inlist,1)-i+1),:) = inlist(i:end,:);




end  % end switch



return