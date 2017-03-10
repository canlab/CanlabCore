function [dt,algo] =  training(algo,dt)

 %% this is a dummy function, because knn is not able to learn any parameter 

   if sum(sum(floor(get_y(dt))==get_y(dt)))~=prod(size(get_y(dt))),
     algo.algorithm.use_signed_output=0;
   end

       disp(['running '  get_name(algo) '.... '])

  algo.dat=dt;
  
  if algo.no_train==1
    dt=set_x(dt,get_y(dt)); 
  else
    dt=test(algo,dt);
  end
  
  
