function w = get_w2(sv,method) 
% w = get_w(sv)
% High value of w means big decrease in the margin.
% Warning: do not return values of w here. Return a score for each feature
% use the max, change max into sum should not change a lot.
% can be used only if same kernel for all base machines
%% method = 1 ->  take the sum
%% method = 2 ->  used for the l0 update, where the get_w is different from rfe
if nargin < 2,
    method = 1;
end;
if method<2,
    
tmp = sv.child{1};
ke=tmp.child;
for i=1:length(sv.child),   
    tmp = sv.child{i};
    if (tmp.child==ke)|strcmp(ke.ker,'weighted_linear'),        
        alphatmp(:,i) = tmp.alpha;
    else
        error('Get w for one-vs-rest valid only when the kernel is the same for all base machines.');
        return;
    end;
end;
Q = size(alphatmp,2);
kertmp = ke.kerparam;
Ktmp = get_kernel(ke,tmp.Xsv,tmp.Xsv);
if strcmp(ke.ker,'linear')|strcmp(ke.ker,'weighted_linear'),
  % compute the multi-class margin for each feature removed.
    if strcmp(ke.ker,'weighted_linear'),
        Xtmp = get_x(tmp.Xsv);
        Xtmp = Xtmp.*(ones(size(Xtmp,1),1)*ke.kerparam);    
        Ktmp = Xtmp*Xtmp';
    end;
    Xtmp = get_x(tmp.Xsv);
    tmp=Xtmp'*alphatmp;
    tmp2=Xtmp'*sum(alphatmp,2);
    w = sum(tmp.^2,2)' + 1/Q*(tmp2.^2)';
    return;
else
    if strcmp(ke.ker,'rbf'),       
      % compute the kernel matrix for all components
      
      Xtmp = get_x(tmp.Xsv);
      
      tmp=zeros(1,Q);
      for j=1:Q,
          tmp(j)=alphatmp(:,j)'*(Ktmp)*alphatmp(:,j);                
      end;
      tmp = sum(tmp) + 1/Q*(sum(alphatmp,2)'*Ktmp*sum(alphatmp,2));
                                    
      % compute the margin when one component is removed
        for i = 1:size(Xtmp,2),
           Ki = Xtmp(:,i)*ones(1,size(Xtmp,1)) - ones(size(Xtmp,1),1)*Xtmp(:,i)';
           Ki = Ki.^2;
           Ki = Ki/(2*kertmp^2); 
           Ki = exp(Ki);
           for j=1:Q,
            tmpi(j)=alphatmp(:,j)'*(Ktmp.*Ki)*alphatmp(:,j);                
           end;
            tmpi = sum(tmpi) + 1/Q*(sum(alphatmp,2)'*(Ktmp.*Ki)*sum(alphatmp,2));                
           w(i) = tmpi;              
        end;
        w=max(w)-w;
       return;    
    elseif strcmp(ke.ker,'poly'),
        Xtmp = get_x(tmp.Xsv);        
        
        tmp=zeros(1,Q);
      for j=1:Q,
          tmp(j)=alphatmp(:,j)'*(Ktmp)*alphatmp(:,j);                
      end;  
      tmp = sum(tmp) + 1/Q*(sum(alphatmp,2)'*Ktmp*sum(alphatmp,2));        
        for i = 1:size(Xtmp,2),
           Ki = Xtmp(:,i)*Xtmp(:,i)';
           K_i = (Ktmp - Ki + 1).^(kertmp);         
           for j=1:Q,
            tmpi(j)=alphatmp(:,j)'*(Ki)*alphatmp(:,j);                
           end;
            tmpi = sum(tmpi) + 1/Q*(sum(alphatmp,2)'*(Ki)*sum(alphatmp,2));                
           w(i) = tmpi;              
        end;
        w=max(w)-w;
       return;
    end;% if strcmp(...,'rbf')
end;% if strcmp(...,'linear')
%%% Otherwise, all classic kernels have been computed. Now use the get_kernel generically
if strcmp(ke.ker,'custom'),
    error('Get w not implemented for CUSTOM kernels.');
    return;
end;
% compute the multi-class margin for each feature removed.
     
   Xtmp = get_x(tmp.Xsv);
   tmp=zeros(1,Q);              
   for j=1:Q,
          tmp(j)=alphatmp(:,j)'*(Ktmp)*alphatmp(:,j);                
   end;
   tmp = sum(tmp) + 1/Q*(sum(alphatmp,2)'*Ktmp*sum(alphatmp,2));        
  
    
   for i = 1:size(Xtmp,2),  
       dtmp = data('tmp',[Xtmp(:,1:(i-1)),Xtmp(:,(i+1):size(Xtmp,2))],[]);
       Ki = get_kernel(ke,dtmp,dtmp);             
       for j=1:Q,
            tmpi(j)=alphatmp(:,j)'*(Ki)*alphatmp(:,j);                
       end;
       tmpi = sum(tmpi) + 1/Q*(sum(alphatmp,2)'*(Ki)*sum(alphatmp,2));                
       w(i) = tmpi;                          
   end;   
   w=max(w)-w;
return;
else, %% if method < 2
    %% method for the multiplicative update
    for i=1:length(sv.child),   
        atmp = sv.child{i};
        W(:,i) = abs(get_w(atmp)');
    end;   
    w = sum(W,2);
    w=w';
end;

