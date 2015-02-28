function [d,a] =  training(a,d),


if a.algorithm.verbosity>0
    disp(['training ' get_name(a) '.... ']);
end;

x=d.X; k=a.k;
%rand('state',2);
l = size(x,1);  


if ~isempty(a.mu)    % -- are cluster centers given? --
  mu=a.mu;
else                  
    if( l<k) 
        warning('@kmeans: k> #data points! setting k=#data points')
        k=l;
    end
  r=randperm(l); mu=x(r(1:k),:);
end
oldy = -ones(1,l); d2 = zeros(k,l); change=1; 

loops=0;

while(change>0)
                     %% -- calc distances from centres --
            
  loops=loops+1; if loops>a.max_loops break; end;
                     
  d1=get_distance(a.child,data(x),data(mu));
  [dummy,y] = min(d1);

  for j=1:k
     if sum(y==j)==0 % if no representatives of cluster j
         f=find(oldy==j); 
         if ~isempty(f) y(f(1))=j; end; % force cluster j to still exist    
     end    
  end   

  
                     %% -- now update means mu --		     
  for i=1:k
    if sum(y==i)>0
      mu(i,:)=mean(x(y==i,:));
    end
  end 

                      %% -- stop if no changes  --

  change = sum(y~=oldy); oldy=y;
  if a.algorithm.verbosity>0 
    disp(['(moved ' num2str(change) ' vectors)']); 
  end;
end 

a.y=y;

a.mu=mu;  
d=test(a,d);



% Old version of code, didn't force every cluster to have at least 
% one trainbing example in it: 
% 
% function [d,a] =  training(a,d),
% 
% 
% if a.algorithm.verbosity>0
%     disp(['training ' get_name(a) '.... ']);
% end;
% 
% x=d.X; k=a.k;
% %rand('state',2);
% l = size(x,1);  
% 
% 
% if ~isempty(a.mu)    % -- are cluster centers given? --
%   mu=a.mu;
% else                  
%   r=randperm(l); mu=x(r(1:k),:);
% end
% oldy = -ones(1,l); d2 = zeros(k,l); change=1; 
% 
% while(change>0)
%                      %% -- calc distances from centres --
%                      
%   d1=get_distance(a.child,data(x),data(mu));
%   [dummy,y] = min(d1);
% 
%                      %% -- now update means mu --		     
%   for i=1:k
%     if sum(y==i)>0
%       mu(i,:)=mean(x(y==i,:));
%     end
%   end 
% 
%                       %% -- stop if no changes  --
% 
%   change = sum(y~=oldy); oldy=y;
%   if a.algorithm.verbosity>0 
%     disp(['(moved ' num2str(change) ' vectors)']); 
%   end;
% end 
% 
% a.mu=mu;  
% d=test(a,d);
% 
% 















	
