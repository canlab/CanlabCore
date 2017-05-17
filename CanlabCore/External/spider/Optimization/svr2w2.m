function [bound,grad] = svr2w2 (sigma,X,Y,ker,p1,scales,slacks,opt,var2)
  
 % if ker=='weighted_linear' p1=1; else p1=2; end    
 %% NOTE: ASSUMES AT MOMENT THAT IT IS A POLYNOMIAL or rbf, if using 'quadprog'-- not compatible w/ svmlight  
 
 % chirag's edits: rbf kernel (non-weighted) implemented;
 % cannot use with 'svmlight' optimizer as of yet..
 
  l = size(X,1);
  N=size(X,2);
  Xkeep=X; 
      
  scale=max(scales); % number of scaling factors 
  if scale>0
    if length(scales)==1
      D=repmat(sigma(1:N),1,l);
    else
      D=repmat(sigma(scales),1,l);      
    end
    X= X .* D';
  end

  slack=max(slacks);
  if slacks==0 
    ridge=ones(l,1)*1e-11; 
  else
    if length(slacks)==1 ridge=ones(l,1)*sigma(scale+1); end; 
    if length(slacks)>1 ridge=sigma(scale+slacks); end; 
  end
  switch ker
      case 'weighted_linear'
          K = X*X';
      case 'weighted_poly'
          K = X*X';
          K=(K+1).^p1;
      case 'rbf'
          %% create rbf kernel...
          [n1,n2]=size(X);
          kerneloption=ones(1,n2)*p1;
          metric = diag(1./p1.^2);
          ps = X*metric*X'; 
          [nps,pps]=size(ps);
          normx = sum(X.^2*metric,2);
          normxsup = sum(X.^2*metric,2);
          ps = -2*ps + repmat(normx,1,pps) + repmat(normxsup',nps,1) ; 
          K = exp(-ps/2);
          %%
  end
  
  K= K + diag(ridge) .^2; %(sum(sigma)/N)*eye(l)*ridge;
  
  H = K .* (Y*Y');
  
  % construct the lower/upper bound vector
  lower = [zeros(l,1)];
  C=Inf;
  upper = [C*ones(l,1)];
  A = Y';  b = 0;   c = [-ones(l,1)];   
  
  switch opt
   case {'loqo'}
    [a,b0] = loqo(sparse(H),c,sparse(Y'),0,lower,upper,[],1,0);
  
   case{'quadprog'}
    options = optimset('LargeScale','off');
    options = optimset(options,'MaxIter',400);
    [a,fval,exit,out,lambda] = quadprog(H,c,[],[],Y',0,lower,upper, ...
					[],options);
    b0=lambda.eqlin(1);
   case{'andre'}
    [a,b0] = quadsolve(H,c,Y',0,10000);
  
   case{'default'}  
    [a,b0]=qp_learn(full(K),Y,1e5);  
  
   case{'svmlight'}
    error('not implemented yet for rbf kernels!');
    %if strcmp(get(al.child,'ker'),'linear') ker=1; p1=1; end;
    %if strcmp(get(al.child,'ker'),'poly') ker=2;
    %p1=get(al.child,'kerparam'); p1=p1{1}; end; 
     %if strcmp(get(al.child,'ker'),'gaussian') ker=3; p1= ...
     %get(al.child,'kerparam'); p1=(2*p1{1}^2); 
     %end;
     ker='linear';p1=1;
     X=full(X);
     [a b0 ind] = svmlight(X,Y,1e5,ridge,0,ker,p1,0);
     alpha=zeros(size(X,1),1); alpha(ind+1)=a;
     a=abs(alpha);
  end
   
  % radius calculation 
  
  if var2==0
    switch opt    
     case {'loqo'}
      b=loqo(sparse(2*K),-diag(K),sparse(ones(1,l)),1,zeros(l,1),[],[],1,0);
     case {'quadprog'}
      b=quadprog(2*K,-diag(K),[],[],ones(1,l),1,zeros(l,1),[],[],options);
     case {'andre'}
       error(['calc radius: not implemented for this optimizer, use  svmlight or quadprog']);
     case {'default'}
      error(['calc radius: not implemented for this optimizer, use  svmlight or quadprog']);
     case {'svmlight'}   %% use andre anyway for r2w2 calculation
      [b] = quadsolve(2*K,-diag(K),ones(1,l),1,10000);
    end
    r=b'*diag(K)-b'*K*b;
    r=abs(r);
  else
    r = mean(diag(K))-mean(mean(K));
  end
  
  w = a'*H*a; %w = sum(a); %
  bound=r*w;
  
  XXtemp=X*X'; 
  if p1>1
      XXtemp=XXtemp+1; 
  end    
  
  grad= zeros(scale+slack,1);
  
  for i=1:scale                       %% calculate scaling factors
    Xc = Xkeep(:,find(i==scales));
    Kc = Xc*Xc';
    switch ker
        case {'weighted_poly', 'weighted_linear'}
            Kgrad= 2 * sigma(i) *  (Kc * p1) .* (XXtemp).^(p1-1);
        case 'rbf'
            x_k=Xc;
            nbx=size(x_k,1);
            deriv=(x_k*ones(1,nbx)-ones(nbx,1)*x_k').^2;
            Kgrad = (deriv.*K)./(p1^2);
    end
    gradw= -(a.*Y)'*Kgrad*(a .* Y);
    if var2==0
      gradr= b'*diag(Kgrad)- b'*Kgrad*b;       
    else
      gradr= mean(diag(Kgrad))-mean(mean(Kgrad));
    end
    grad(i)=  (gradw*r + gradr*w);
  end; 
  

  
  
  for i=1:slack                       %% calculate slacking factors
    rs=ridge; rs(find(i~=slacks))=0;
    Kgrad=2*diag(rs);
    gradw= -(a.*Y)'*Kgrad*(a .* Y);
    if var2==0
      gradr= b'*diag(Kgrad)- b'*Kgrad*b;       
    else
      gradr= mean(diag(Kgrad))-mean(mean(Kgrad));
    end
    grad(i+scale)=  (gradw*r + gradr*w);
  end
  






