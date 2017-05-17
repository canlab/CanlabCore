function [retdat,algo] =  training(algo,dat)

retdat=dat;

 disp(['training ' get_name(algo) '.... '])
 switch algo.optimizer
  case {'svmtorch'}
  %%<<----------------svmtorch optimizer-------------------->>

    multi = 0;
    regression = 1;   
    degree = 1;
    gamma = 1; 
    eps = algo.epsilon;
    C = algo.C;
    if strcmp(algo.child.ker,'linear')
       kernelType = 0;
    elseif strcmp(algo.child.ker,'poly')
       kernelType = 1;
       degree = algo.child.kerparam;
    elseif strcmp(algo.child.ker,'rbf')
    	kernelType = 2;    
	gamma = algo.child.kerparam;
    end;
    [x y] = get_xy(dat);
    numEx = get_dim(dat);    
    [alpha,threshold0,xSV] = SVMTorch(x,y,regression,multi,kernelType,degree,gamma,C,eps);
    if isempty(alpha)
    	alpha = zeros(numEx,1);
    end;

 %<<--------------sparse svr------------------>>%

   case{'sparse'}
   	kern1=get_kernel(algo.child,dat,dat);
	y=get_y(dat);
	kNum=size(kern1,1);
	c=ones(kNum*4,1);
	c(2*kNum+1:4*kNum)=algo.C;
	kern2 = [-kern1 kern1 -1*eye(kNum) zeros(kNum); kern1 -kern1 zeros(kNum) -1*eye(kNum); -1*eye(4*kNum)];
	cst = zeros(6*kNum,1);
	cst(1:kNum) = (algo.epsilon+algo.b0) - y;
	cst(kNum+1:2*kNum) = (algo.epsilon-algo.b0) + y;	
	opts= optimset('display','off','MaxIter',10000,'LargeScale','off');
	[alphaTemp,fval,exit,out,lambda] = linprog(c,kern2,cst,[],[],[],[],[],opts);
	threshold0=lambda.ineqlin(1);
	alpha=alphaTemp(1:kNum)-alphaTemp(kNum+1:2*kNum);

  %<<------------andre optimizer--------------------->> 

  case {'andre'}
   kern1=get_kernel(algo.child,dat,dat);   %% <--- calculate kernel
   y=get_y(dat); 
   yLen = length(y(:,1));
   xSV = get_x(dat);
   kern2 = [kern1 , -kern1 ; -kern1 , kern1];   
   cst = ones(2*yLen,1); 
   cst(yLen+1:2*yLen) = -1;
   if algo.nu ==0,
       c = zeros(2*yLen,1);
       c(1:yLen) = algo.epsilon*ones(yLen,1) - y;
       c(yLen+1:2*yLen) = algo.epsilon*ones(yLen,1) + y;       
       [alphaTemp,threshold] = quadsolve(kern2,c,cst',0,algo.C); 
       alpha=alphaTemp(1:yLen)-alphaTemp(yLen+1:2*yLen);
       threshold0 = -threshold;
   else
     if algo.C==Inf,
         algo.C=10000;
     end;

        c = zeros(2*yLen+1,1);
        c(1:yLen) = -y;
        c(yLen+1:2*yLen) = y;
        cst2 = ones(1,2*yLen+1)/(yLen*algo.nu);
        cst2(2*yLen+1)=-1;
        kern2 = [kern2,zeros(2*yLen,1);zeros(1,2*yLen+1)];
        cst=[cst',0;cst2];
        [alphaTemp,threshold] = quadsolve(kern2,c,cst,[0;0],algo.C); threshold0 = -threshold(1);
        alpha=alphaTemp(1:yLen)-alphaTemp(yLen+1:2*yLen);
        epsilon = -threshold(2);
    end
    
  %<<------------quadprog optimizer--------------------->> 
  
  case {'quadprog'}
   
   kern1=get_kernel(algo.child,dat,[]);   %% calculate kernel
   kNum=size(kern1,1);  y=get_y(dat); 
   kern2 = [kern1 , -kern1 ; -kern1 , kern1];
   c = zeros(2*kNum,1);
   c(1:kNum) = algo.epsilon*ones(kNum,1) - y;
   c(kNum+1:2*kNum) = algo.epsilon*ones(kNum,1) + y;
   cst(kNum+1:2*kNum) = -1;
   cst = ones(2*kNum,1); 
   opts= optimset('display','off','MaxIter',10000,'LargeScale','off'); 
   [alphaTemp,fval,exit,out,lambda] = quadprog(kern2,c,[],[],cst',0,...
				       zeros(2*kNum,1),algo.C*ones(2*kNum,1),[],opts);
   threshold0=lambda.eqlin(1);
   alpha=alphaTemp(1:kNum)-alphaTemp(kNum+1:2*kNum);
  
%<<------------libsvm optimizer--------------------->> 
case {'libsvm'}
  
        %      
        x=[];
        y=[];
        svm_type=3;
        kernelType=0;
        degree=3;
        gamma=0;
        coef0=0;
        
        nu=algo.nu;
        if(nu>0)
            svm_type=4;
        end

        
        cachesize=40;
        C=algo.C;
        eps=algo.epsilon;
        p=0.05;
        shrinking=1;
        
        
        weight_label=[];
        weight=[];
        nr_weight=0;
        
        if strcmp(algo.child.ker,'linear')
            kernelType = 0;
        end;
        if strcmp(algo.child.ker,'poly')
            kernelType = 1; 
            degree = algo.child.kerparam;
            coef0 = 1;
            gamma = 1;
        end;
        if strcmp(algo.child.ker,'rbf'),
            kernelType = 2; 
            sigma = algo.child.kerparam; 
            gamma = 1/(2*sigma^2);
        end;

        y=get_y(dat); 
        x=get_x(dat);

        if strcmp(algo.child.ker,'custom'),
          kernelType = 4; 
          K= algo.child.kerparam;
          l = get_dim( retDat);
          x = get_index( retDat);
          x = [ reshape( x, l, 1) [ 1:l]']; % using x to pass indices in Matrix and real indices
        end;
 
        s=whos('libsvm_cachesize','global');
        
        if (length(s)>0)
            global libsvm_cachesize;
            cachesize=libsvm_cachesize;
        else
            cachesize=40;
        end
        if algo.algorithm.verbosity>1
         fprintf('Using %d MB Cache for Libsvm\n',cachesize)
        end

    
         if( kernelType == 4)
          [alpha,xSV,bias0]=libsvm_regressor_spider(x,y,svm_type,kernelType,...
                     degree,gamma,coef0,nu,cachesize,C,eps,p,weight_label,weight,nr_weight,K);
          algo.Xsv=get(retDat, xSV( :, 2));
        else
          [alpha,xSV,bias0]=libsvm_regressor_spider(x,y,svm_type,kernelType,...
                     degree,gamma,coef0,nu,cachesize,C,eps,p,weight_label,weight,nr_weight);
          algo.Xsv=data(xSV);
        end
       

        
        threshold0 = bias0 * y(1); 
        
        algo.b0=bias0;

        algo.epsilon = eps;
        
        algo.Xsv = data(xSV);
        algo.alpha=alpha;
      

        if algo.algorithm.do_not_evaluate_training_error==1   
            retdat=set_x(dat,get_y(dat));
        else
            retdat=test(algo,dat);
        end
        
        return

        
        %         fin=find(abs(alpha)>algo.alpha_cutoff);

% case {'libsvm'}
%   
%    if algo.nu ==0,
%         svm_type = 3; 
%         C = algo.C; 
%         epsilon = algo.epsilon;  
%         nu=0;
%     else
%         svm_type = 4; 
%         nu = algo.nu; 
%         C = algo.C; 
%         epsilon = -1;
%     end;
%     %% default values for libsvm
%     cacheSize = 40; 
%     eps = 0.001; 
%     shrinking=1;
%     nrWeight = 0; 
%     weightLabel =0; 
%     weight = 1; 
%     gamma=1; 
%     deg = 0; 
%     coef0 = 0; 
%     kerTmp = algo.child;
%     if strcmp(kerTmp.ker,'linear')
%       kernelType = 0;
%     end;
%     if strcmp(kerTmp.ker,'poly')
%       kernelType = 1; 
%       ptmp = kerTmp.kerparam; 
%       deg = ptmp; 
%       coef0 = 1;
%     end;
%    if strcmp(kerTmp.ker,'rbf'),
%          kernelType = 2; 
%          ptmp = kerTmp.kerparam; 
%          gamma = 1/(2*ptmp^2);
%    end;
% %    if algo.balanced_ridge~=0,
% %         disp('Warning: balanced ridge not implemented for libsvm.');
% %    end;
%    y=get_y(dat); 
%    x=get_x(dat);
%    [alpha,threshold0,xSV,eps,CC] = svmlibtrain(x,y,svm_type,kernelType,deg,gamma,coef0,nu,cacheSize,C,eps,epsilon,...
%        shrinking,nrWeight,weightLabel,weight,0);    
%    threshold0=-threshold0;% from libsvm
%    % alpha is reordered in order to have the same xsp for all runs (important for one_vs_rest)
%    alphaTemp = zeros(size(xSV,1),1);
%    indTemp = find(xSV(:,size(xSV,2))~=0);
%    indTemp2 = xSV(indTemp,size(xSV,2));
%    alphaTemp(indTemp2) = alpha(indTemp);
%    alpha = alphaTemp;
%    epsilon=eps;
%  end
end

 if algo.nu~=0,
   algo.epsilon = epsilon;
 end;

 
 algo.b0=threshold0;
 fin=find(abs(alpha)>algo.alpha_cutoff);
 algo.alpha=alpha(fin);
 algo.Xsv = get(dat,fin);
      


 if algo.algorithm.do_not_evaluate_training_error==1   
   retdat=set_x(dat,get_y(dat));
 else
   retdat=test(algo,dat);
 end
 
