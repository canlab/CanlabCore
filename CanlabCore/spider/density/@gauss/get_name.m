function s=get_name(a)

s=[get_name(a.algorithm)];  %% print name of algorithm

if a.l~=50
  s=[s ' l=' num2str(a.l) ];             
end
if ~isempty(a.assume)
  s=[s ' assume=' a.assume ];             
end


eval_name                   %% print values of hyperparameters
