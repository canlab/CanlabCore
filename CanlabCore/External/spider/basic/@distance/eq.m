function out=eq(k1,k2),
% check equality between kernels
out=0;
if isempty(k1)&isempty(k2),
   out=1;
return;
end;
if isempty(k1)|isempty(k2),
  out=0;
return;
end;
if ~strcmp(k1.ker,k2.ker),
  return;
else
if iscell(k1.kerparam)&iscell(k2.kerparam)
     if length(k1.kerparam)==length(k2.kerparam),
       for i=1:length(k1.kerparam),	      
               tmp1 = size(k1.kerparam{i});
               tmp2 = size(k2.kerparam{i});
               if tmp1(1)~=tmp2(1)|tmp1(2)~=tmp2(2),
		        return;
               end
               if ~isempty(k1.kerparam{i})|~isempty(k2.kerparam{i}),                   
	            t = sum(sum(k1.kerparam{i}==k2.kerparam{i}));
                if t~=prod(tmp1),
		         return;
	            end 
               end
       end       
       out=1;
     end
else
     if ~iscell(k1.kerparam)&~iscell(k2.kerparam),
        tmp1 = size(k1.kerparam);
        tmp2 = size(k2.kerparam);
        if tmp1(1)~=tmp2(1)|tmp1(2)~=tmp2(2),return;end;
        if isempty(k1.kerparam)&isempty(k2.kerparam), out=1; return; end;
	    t = sum(sum(k1.kerparam==k2.kerparam));
        if t~=prod(tmp1),return;end;       
        out=1;
     end
end
end;
