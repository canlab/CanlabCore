function dt =  testing(algo,dt)
% 
%   Generates labels using knn.
%    
  [numOfSam,vecDim,outDim]=get_dim(dt);


  % kernel or distance - that is the question ?!
  if isa( algo.child, 'distance')
    distflag = 1;
  else
    distflag = 0;
    trainXP=get_norm(algo.child,algo.dat).^2;
    testXP=get_norm(algo.child,dt).^2;
  end

  
%  flag=0;
%  if isa(algo.child,'distance')
%    if strcmp(algo.child.dist,'custom') 
%      flag=1;
%    end
%    if strcmp(algo.child.dist,'custom_fast') 
%      flag=1;
%    end
%  end
%  if flag==0 %% ok to compute norm 
%   trainXP=get_norm(algo.child,algo.dat).^2;
%   testXP=get_norm(algo.child,dt).^2;
%  end
  
    
  y=get_y(algo.dat);
  yTemp=get_y(dt);
  yEst=yTemp*0-1;

% Possible to calc in batch or as for loop, whereas the batch way requires
% the whole kernel matrix to fit in memory.
if algo.batch==1,
  
  if distflag==0
    dist = trainXP*ones(1,length(testXP)) - 2*calc(algo.child,dt,algo.dat) + ones(length(trainXP),1)*testXP';
  else
    dist = calc(algo.child,dt,algo.dat);
  end
  
  [val,valInd]=sort(dist); 
 
  if(~isempty(y))
      for i=1:numOfSam,
          yTmp(i,:) = mean(y(valInd(1:algo.k,i),:),1); 
      end;
  else
   yTmp=[];    
  end        
 if algo.algorithm.use_signed_output==1,
   if size(y,2)>1  %% <---multi-class pattern recognition
     [maxComp maxCompInd]=max(yTmp,[],2); 
     for i = 1:numOfSam,
       yEst(i,maxCompInd(i))=1;
     end;
   else    %%<--- binary pattern recognition
     yEst=sign(yTmp);
   end
 else
   yEst=yTmp; 
 end
 
 if algo.output_preimage 
   yEst=[]; ind=algo.dat.index;
   for i=1:algo.k
     yEst=[yEst ind(valInd(i,:))'];
   end
 end 
 
else
%<<-----online calculations, instead of batch (does not need as much memory
%           as batch)----->
    
  for i=1:numOfSam 
    dist= (trainXP - 2*get_kernel(algo.child,dt,algo.dat,i) ) + testXP(i);
    
    [val,valInd]=sort(dist); 
   
   if algo.algorithm.use_signed_output==1,
    if size(y,2)>1  %% <---multi-class pattern recognition 
      tmp=sum(y(valInd(1:algo.k),:),1); 
      [maxComp maxCompInd]=max(tmp); 
      yEst(i,maxCompInd)=1;
    else            %% binary pattern recognition
      yEst(i)=sign(mean(y(valInd(1:algo.k),:)));
    end
   else
     yEst(i,:)=mean(y(valInd(1:algo.k),:),1); 
   end 
  end
  
  if algo.output_preimage 
    ind=algo.dat.index;
    for j=1:algo.k
      yEst(i,j)=[yEst ind(valInd(j,:))'];
    end
   end;
  
%<<-----------end online------------------>>
end;
  
  if algo.algorithm.use_signed_output==1 %% <--- break deadlocks
    yEst(yEst==0)=-1;
  end

  dt=set_x(dt,yEst);
  dt=set_name(dt,[get_name(dt) ' -> ' get_name(algo)]); 

