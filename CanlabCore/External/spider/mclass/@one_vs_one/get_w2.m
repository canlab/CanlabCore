function weig = get_w2(sv,method) 
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
    

    temp = sv.child{1};
    ke=temp.child;
    for i=1:length(sv.child),   
        temp = sv.child{i};
            if (temp.child==ke)|strcmp(ke.ker,'weighted_linear'),        
              alphaTemp(:,i) = temp.alpha;
            else
                error('Get w for one-vs-rest valid only when the kernel is the same for all base machines.');
            return;
            end;
    end;
    Q = size(alphaTemp,2);
    kertmp = ke.kerparam;
    Ktmp = get_kernel(ke,temp.Xsv,temp.Xsv);
    if strcmp(ke.ker,'linear')|strcmp(ke.ker,'weighted_linear'),
  % compute the multi-class margin for each feature removed.
        if strcmp(ke.ker,'weighted_linear'),
         xTemp = get_x(temp.Xsv);
         xTemp = xTemp.*(ones(size(xTemp,1),1)*ke.kerparam);    
         Ktmp = xTemp*xTemp';
        end;
        xTemp = get_x(temp.Xsv);
        temp=xTemp'*alphaTemp;
        tmp2=xTemp'*sum(alphaTemp,2);
        weig = sum(temp.^2,2)' + 1/Q*(tmp2.^2)';
        weig = max(weig)-weig;
        return;
    else
        if strcmp(ke.ker,'rbf'),       
      % compute the kernel matrix for all components
      
         xTemp = get_x(temp.Xsv);
      
        temp=zeros(1,Q);
        for j=1:Q,
          temp(j)=alphaTemp(:,j)'*(Ktmp)*alphaTemp(:,j);                
        end;
        temp = sum(temp) + 1/Q*(sum(alphaTemp,2)'*Ktmp*sum(alphaTemp,2));
                                    
      % compute the margin when one component is removed
        for i = 1:size(xTemp,2),
           Ki = xTemp(:,i)*ones(1,size(xTemp,1)) - ones(size(xTemp,1),1)*xTemp(:,i)';
           Ki = Ki.^2;
           Ki = Ki/(2*kertmp^2); 
           Ki = exp(Ki);
           for j=1:Q,
            tmpi(j)=alphaTemp(:,j)'*(Ktmp.*Ki)*alphaTemp(:,j);                
           end;
            tmpi = sum(tmpi) + 1/Q*(sum(alphaTemp,2)'*(Ktmp.*Ki)*sum(alphaTemp,2));                
           weig(i) = temp-tmpi;              
        end;
        weig = max(weig)-weig;
       return;    
    elseif strcmp(ke.ker,'poly'),
        xTemp = get_x(temp.Xsv);        
        
        temp=zeros(1,Q);
      for j=1:Q,
          temp(j)=alphaTemp(:,j)'*(Ktmp)*alphaTemp(:,j);                
      end;
      temp = sum(temp) + 1/Q*(sum(alphaTemp,2)'*Ktmp*sum(alphaTemp,2));        
        for i = 1:size(xTemp,2),
           Ki = xTemp(:,i)*xTemp(:,i)';
           K_i = (Ktmp - Ki + 1).^(kertmp);         
           for j=1:Q,
            tmpi(j)=alphaTemp(:,j)'*(Ki)*alphaTemp(:,j);                
           end;
            tmpi = sum(tmpi) + 1/Q*(sum(alphaTemp,2)'*(Ki)*sum(alphaTemp,2));                
           weig(i) = temp-tmpi;              
        end;
       weig=max(weig)-weig;
       return;
    end;% if strcmp(...,'rbf')
end;% if strcmp(...,'linear')
%%% Otherwise, all classic kernels have been computed. Now use the get_kernel generically
if strcmp(ke.ker,'custom'),
    error('Get w not implemented for CUSTOM kernels.');
    return;
end;
% compute the multi-class margin for each feature removed.
     
   xTemp = get_x(temp.Xsv);
   temp=zeros(1,Q);              
   for j=1:Q,
          temp(j)=alphaTemp(:,j)'*(Ktmp)*alphaTemp(:,j);                
   end;
   temp = sum(temp) + 1/Q*(sum(alphaTemp,2)'*Ktmp*sum(alphaTemp,2));        
  
    
   for i = 1:size(xTemp,2),  
       datTemp = data('temp',[xTemp(:,1:(i-1)),xTemp(:,(i+1):size(xTemp,2))],[]);
       Ki = get_kernel(ke,datTemp,datTemp);             
       for j=1:Q,
            tmpi(j)=alphaTemp(:,j)'*(Ki)*alphaTemp(:,j);                
       end;
       tmpi = sum(tmpi) + 1/Q*(sum(alphaTemp,2)'*(Ki)*sum(alphaTemp,2));                
       weig(i) = temp-tmpi;                          
   end;   
   weig=max(weig)-weig;
return;
else, %% if method < 3
    %% method for the multiplicative update
    for i=1:length(sv.child),   
        atmp = sv.child{i};
        W(:,i) = abs(get_w(atmp)');
    end;   
    weig = sum(W,2);
    weig=weig';
end;
