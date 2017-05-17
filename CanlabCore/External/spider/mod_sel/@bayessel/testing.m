function d =  testing(a,d)  

  d=set_x(d,a.pbest); 
  d=set_y(d,[]);
  %s = ['gamma = ' num2str(a.pbest(1)) ' C = ' num2str(a.pbest(2))];
  %if length(a.pbest) == 3 % regression
%	s = [s ' nu = ' num2str(a.pbest(3))];
 % end;
  s = ['C = ' num2str(a.pbest(1))];
  d=set_name(d,[get_name(a) ' --> ' s]); 
