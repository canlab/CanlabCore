function [max,maxrest] = getMaxInARow(stimList,varargin)
	if nargin>1,dojitter = varargin{1},else dojitter = 0;,end

	max = 1; counter = 1;maxrest = 1;counterrest = 1;
        if dojitter, restnum = max(stimList);,else restnum = Inf;,end
	for i = 2:size(stimList,1)
		if stimList(i,1) == stimList(i-1,1)& ~stimList(i,1) == 0,
			counter = counter + 1;
			if counter > max, max = counter;,end
		else counter = 1;
                end
                
                if dojitter
                     if stimList(i,1) == restnum & stimList(i-1,1) == restnum
                        counterrest = counterrest+1;
                        if counterrest > maxrest,maxrest = counterrest;,end
                     else counterrest = 1;
		     end
                end

	end

return


