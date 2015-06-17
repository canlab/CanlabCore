function sl = transform2block(sl,rl,restlength,dofirst)
% sl and rl are stimlist and restlist.
% now is NOT compatible with a separate jitter blank rest interval.
% insert rests before transforming with this function.
%
% dofirst: transforms 1st trial of each block into its OWN trial type.
% doesn't work right now - transforms only every other block.

if isempty(restlength), restlength = 0;,end
origlength = length(sl);

maxnum = max(sl);
%if doblank, maxnum = maxnum-1;,end

for i = 2:2:length(rl)
     start = sum(rl(1:i-1)) + restlength * length(rl(1:i-1)) + 1;
     stop = min(start + rl(i) - 1,length(sl));
     sl(start:stop) = (sl(start:stop) > 0) .* (sl(start:stop) + maxnum);
     %if dofirst
	%	sl(start) = 2 * maxnum + 1;
    %end
end

if dofirst
    for i = 1:length(rl)
        start = sum(rl(1:i-1)) + restlength * length(rl(1:i-1)) + 1;
        sl(start) = 2 * maxnum + 1;
    end
end
        
sl = sl(1:origlength);

return
