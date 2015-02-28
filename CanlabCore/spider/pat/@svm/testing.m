function dat  =  testing(algo,dat)

%  if isempty(algo.child.dat) % maybe used a custom ker to train and then changed!
%    kerMaTemp=get_kernel(algo.child,dat,algo.xSV);
%  else
%    kerMaTemp=test(algo.child,dat);  
%    kerMaTemp=kerMaTemp.X;
%  end

  %% <<---To avoid large kernel matrices, we test in batches---> 
  sz=get_dim(algo.Xsv);   %% <---500x500 point are the maximum for one batch
  if sz==0 sz=1; end;
  batch_size=round((500^2)/sz);
   
  yEst=[];
  for i=[1:batch_size:get_dim(dat)]
    take= [i:min(i+batch_size-1,get_dim(dat))];
    if ~isempty(algo.alpha) 
     kerMaTemp=get_kernel(algo.child,get(dat,take),algo.Xsv);
     yEst=[yEst; ((algo.alpha'* kerMaTemp)+algo.b0)'];
    else
        yEst=[yEst ; algo.b0*ones(length(take),1)];    
    end
  end
  
  if algo.algorithm.use_signed_output==1
    yEst=sign(yEst);
  end

  dat=set_x(dat,yEst); 
  dat=set_name(dat,[get_name(dat) ' -> ' get_name(algo)]); 
 