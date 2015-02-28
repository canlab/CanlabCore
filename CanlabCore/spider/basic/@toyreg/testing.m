function [d] =  testing(a,d)
  
  %name=[]; if ~isempty(d) name=d.name; end; % record old name
  d=generate(a);  
  %d.name=[name ' -> ' d.name];
  

  