function dat = calc(e,dat)
   if isa(dat,'cell')
     for i=1:length(dat)
       dat{i}=feval(e.type,e,dat{i});
     end
   else
     dat=feval(e.type,e,dat);
   end
   