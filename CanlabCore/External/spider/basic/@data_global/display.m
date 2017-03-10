function [s,it]=display(c,tab,it)
  
  if nargin>1
    [s,it]=display(c.algorithm,tab,it);
    return;
  end
  
  disp(c.algorithm.name);
  disp(['data dimensions:'])
  
  global X; global Y;
   
  if isempty(c.myX)
    s=(['X = ' num2str(size(X(c.index,c.findex),1)) 'x' ...
	  num2str(size(X(c.index,c.findex),2))  ]);
  else 
    s=(['X = ' num2str(size(c.myX,1)) 'x' num2str(size(c.myX,2))  ]);
  end
  
  
  if isempty(c.myY) 
    if ~isempty(Y)
      disp([s '  Y = ' num2str(size(Y(c.index,:),1)) 'x' ...
	    num2str(size(Y(c.index,:),2))  ]);
    else
      disp([s '  Y = 0x0']);
    end
  else
    disp([s '  Y = ' num2str(size(c.myY,1)) 'x' num2str(size(c.myY,2))  ]);
  end
