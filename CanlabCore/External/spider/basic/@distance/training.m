
function [d,a] = training(a,d)

  a.dat=d;
  d=set_x(d,feval('get_distance',a,d)); 

