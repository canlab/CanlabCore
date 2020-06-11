function [reslt,tReslt,alg] =  traintests(algo,dat,trn,tst,lossType)
  

% <<-----unspecified loss ------>>
  if 
      nargin==4 lossType=[];  
  end  
      
  if ~isa(dat,'cell')
    [dat1,dat2,algo]=traintest(algo,dat,trn,tst,lossType);
  else     
    for i=1:length(dat)
      [res1,res2,algo]=traintest(algo,dat{i},trn,tst,lossType);
      dat1{i}=res1; 
      dat2{i}=res2;
    end
  end
  reslt=dat1;
  tReslt=dat2;
  alg=algo;
