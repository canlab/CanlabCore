function res =  testing(alg,dat,lossType)  
  
  dat.name = [dat.name ' -> ' alg.algorithm.name ];  
  res=test(alg.best,dat);  
  
  
