function ret=eq(k1,k2),
%
% check equality between kernels
%
ret=0;
if isempty(k1)&isempty(k2),
   ret=1;
return;
end;
if isempty(k1)|isempty(k2),
  ret=0;
return;
end;
if ~strcmp(k1.ker,k2.ker),
  return;
else
if iscell(k1.kerparam)&iscell(k2.kerparam)
     if length(k1.kerparam)==length(k2.kerparam),
       for i=1:length(k1.kerparam),	      
               temp1 = size(k1.kerparam{i});
               temp2 = size(k2.kerparam{i});
               if temp1(1)~=temp2(1)|temp1(2)~=temp2(2),
		        return;
               end
               if ~isempty(k1.kerparam{i})|~isempty(k2.kerparam{i}),                   
	            sm = sum(sum(k1.kerparam{i}==k2.kerparam{i}));
                if sm~=prod(temp1),
		         return;
	            end 
               end
       end       
       ret=1;
     end
else
     if ~iscell(k1.kerparam)&~iscell(k2.kerparam),
        temp1 = size(k1.kerparam);
        temp2 = size(k2.kerparam);
        if temp1(1)~=temp2(1)|temp1(2)~=temp2(2),
            return;
        end;
        if isempty(k1.kerparam)&isempty(k2.kerparam), 
            ret=1; 
            return; 
        end;
	    sm = sum(sum(k1.kerparam==k2.kerparam));
        if sm~=prod(temp1),
            return;
        end;       
        ret=1;
     end
end
end;
