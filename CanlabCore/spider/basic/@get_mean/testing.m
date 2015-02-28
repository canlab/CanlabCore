
function r = testing(a,d)
   
    
  if isa(a.child,'algorithm')
    d=test(a.child,d);
  end
  
  d.group='separate'; %% temporarily deal with each item separately
  if isa(d,'data') r=d; return; end;
  
  r=calc_mean(d,a.loss_type,0,a.take_average);
    
