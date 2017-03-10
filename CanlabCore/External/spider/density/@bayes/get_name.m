function s=get_name(a)

s=[get_name(a.algorithm)];  %% print name of algorithm

if iscell(a.child)
  s=[s '(' get_name(a.child{1}) ',...)'];             
else
  s=[s '(' get_name(a.child) ')'];             
end

eval_name                   %% print values of hyperparameters
