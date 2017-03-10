function [results,a] =  training(a,d)
       
  % [results,algorithm] =  training(algorithm,data,loss)
  disp(['training ' get_name(a) '.... '])

  %% Perform a cv
  s = ['folds=' num2str(a.folds)];
  a.child.use_signed_output=0;
  
  algtmp = cv(a.child,s);
  r = train(algtmp,d);
  
  %% Get the results
  tmpX=[];
  tmpY=[];
  for i=1:length(r.child),
      tmpX =[tmpX;r{i}.X];       
      tmpY = [tmpY;r{i}.Y];
  end;
  
  %% Find the coeff of the sigmoid
  [A,B] = platt_opt(tmpX,tmpY,length(find(tmpY==-1)),length(find(tmpY==1)));  
  a.A=A;
  a.B=B;  

  %% Retrain the child on the whole training set
  [r,a.child] = train(a.child,d);

  results=test(a,d);


 %%% Auxiliary functions


  function [A,B]=platt_opt(out,target,prior0,prior1)
% find the coefficients A and B such that the posterior probability
% of P(y=1|x) = 1/(1+exp(A*f(x)+B)), where f(x) is the output
% of the SVM
%
% out: vector of outputs of the SVM on a validation set
% target: validation labels
% prior0: number of negative points
% prior1: number of positive points
%
% If no validation set is available, one might use the training
% set to do a leave-one-out procedure. Using the span, this means
% replacing out by something like out-target.*alpha.*span

  A=0;
  B=log((prior0+1)/(prior1+1));
  hiTarget= (prior1+1)/(prior1+2);
  loTarget= 1/(prior0+2);
  lambda=1e-3;
  olderr=1e300;
  pp=ones(length(out),1)*(prior1+1)/(prior0+prior1+2);
  count=0;
  t = (target==1)*hiTarget + (target==-1)*loTarget;
  
  for it=1:100
    d1=pp-t;
    d2=pp.*(1-pp);
    a=sum(out.*out.*d2);
    b=sum(d2);
    c=sum(out.*d2);
    d=sum(out.*d1);
    e=sum(d1);
    if (abs(d) < 1e-9 & abs(e) < 1e-9)
      break;
    end;
    oldA=A;
    oldB=B;
    err=0;
    while (1)
      deter=(a+lambda)*(b+lambda)-c*c;
      if deter==0 
	lambda=lambda*10;
	break;    
      end;
      A=oldA+((b+lambda)*d-c*e)/deter;
      B=oldB+((a+lambda)*e-c*d)/deter;
      
      pp=(1+exp(out*A+B)).^(-1);
      pp2=(1+exp(-out*A-B)).^(-1);
      err = -sum(t.*log(pp)+(1-t).*log(pp2));  
      
      if err < olderr*(1+1e-7)
	lambda=lambda*0.1;
	break;
      end;
      
      
      lambda=lambda*10;
      if lambda>1e6
	error('Lambda too big');
      end;
    end;
    
    diff=err-olderr;
    scale=0.5*(err+olderr+1);
    if (diff > -1e-3*scale) & (diff < 1e-7*scale)
      count=count+1;
    else
      count=0;
    end;
    olderr=err;
    if count==3
      break;
    end;
  end;






