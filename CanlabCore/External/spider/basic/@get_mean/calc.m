function d = calc(e,d)

if 0
   if isa(d,'cell')
     for i=1:length(d)
       d{i}=feval(e.loss_type,e,d{i});
     end
   else
     d=feval(e.loss_type,e,d);
   end
end