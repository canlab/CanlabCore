function [ans]= myexist(which,type)

ans=0;

if 0   %% this caused a huge slow down on some machines, and 
       %% is rarely used anyway  
  ans=exist(which,type);
  if strcmp(which,'ridge') ans=0; end;
  if strcmp(which,'poly') ans=0; end;
  if strcmp(which,'alpha') ans=0; end;
  if strcmp(which,'rank') ans=0; end;
  if strcmp(which,'sigma') ans=0; end;
end
