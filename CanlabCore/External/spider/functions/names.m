
function [n,r] = names(d,loss_type)
  
  if nargin==1 loss_type=[]; end;

  if isa(d,'group')
    d=group2cell(d);
  end
  
  r=[];
  d=make_cell(d);
  for i=1:length(d)
    n{i}=[num2str(i) ':' get_name(d{i})];
    
    if ~isempty(loss_type)
      l=loss(d{i},loss_type);
      if prod(size(l))==1
	n{i}=[n{i}  ' -> ' loss_type ' = ' num2str(l) ];
      end
      r=[r l];
    else
      if prod(size(get_y(d{i})))==1 %% make array of losses 
	r=[r get_y(d{i})];
      end
    end
  end
  n=char(n);
