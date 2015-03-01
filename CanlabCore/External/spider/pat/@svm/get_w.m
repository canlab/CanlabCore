function weights = get_w(supVec,method) 

%
%               w = get_w() 
% 
%   Returns the weight vector from the svm 
%   AE: i have changed this procedure to include non linear feature
%   selection with RFE. When non-linear, it gives  a score vector for each
%   feature (used in RFE) which is not the margin but is sorted as the margin


kerTemp = supVec.child.kerparam;
if strcmp(supVec.child.ker,'linear')
  weights=supVec.alpha'*get_x(supVec.Xsv);
  return;
end
 
if strcmp(supVec.child.ker,'weighted_linear')
  len=get_dim(supVec.Xsv);
  weights=supVec.alpha'* (get_x(supVec.Xsv) .*  repmat(kerTemp,len,1));
return;
end
if strcmp(supVec.child.ker,'rbf')  | strcmp(supVec.child.ker,'poly'),
    xTemp = get_x(supVec.Xsv);
    alphaTemp = supVec.alpha;
    
    if strcmp(supVec.child.ker,'rbf'),
        
      % compute the kernel matrix for all components
      kernelMatrix = xTemp*xTemp';
      kerMaDegN = sum(xTemp.^2,2);
      kernelMatrix = ones(size(xTemp,1),1)*kerMaDegN' + kerMaDegN*ones(1,size(xTemp,1)) - 2*kernelMatrix;
      kernelMatrix = exp(-kernelMatrix/(2*kerTemp));
      % compute the margin when one component is removed
       for i = 1:size(xTemp,2),
           kerMai = xTemp(:,i)*ones(1,size(xTemp,1)) - ones(size(xTemp,1),1)*xTemp(:,i)';
           kerMai = kerMai.^2;
           kerMai = kerMai/(2*kerTemp); 
           kerMai = exp(kerMai);
           weights(i) = (alphaTemp'*(kernelMatrix.*kerMai)*alphaTemp);
       end;
    elseif strcmp(supVec.child.ker,'poly'),
        % compute the margin when one component is removed
        Ktmp = xTemp*xTemp';
        for i = 1:size(xTemp,2),
           kerMai = xTemp(:,i)*xTemp(:,i)';
           K_i = (Ktmp - kerMai+1).^(kerTemp);           
           weights(i) = (alphaTemp'*K_i*alphaTemp);
    	end;    
    end;% if strcmp(...,'rbf')
    weights=max(weights)-weights;  %% make largest the smallest --- wrong way round!
    return;
end;% if strcmp(... | strcmp(...,'poly')
%%% if this point is reached then the kernel is a non classic one
ker=supVec.child.ker;
if strcmp(ker,'custom'),
    error('Get w not implemented for CUSTOM kernels.');
    return;
end;
%%% if tpoly kernel -> used for nfe
if strcmp(ker,'tpoly'),
    [len,vecDim,outDim]=get_dim(supVec.Xsv);    
    xTemp = get_x(supVec.Xsv);
    alphaTemp = supVec.alpha;
    weights = ones(length(supVec.child.param,vecDim));
    
    for j=1:length(supVec.child.kerparam),  
        ktmp = supVec.child;
        kerp = ktmp.kerparam;
        for i = 1:length(kerp{j}),
            kerptmp=kerp;
            kerptmp{j} = [kerp{j}(1:i-1),kerp{j}(i+1:length(kerp{j}))];    
            ktmp.kerparam=kerptmp;
            K_i = get_kernel(ktmp,supVec.Xsv,supVec.Xsv);            
            weights(j,kerp{j}(i)) = (alphaTemp'*K_i*alphaTemp);
        end; 
    end;
    weights = max(max(weights))-weights;
    return;
end;

xTemp=supVec.Xsv;  alphaTemp=supVec.alpha; 
[l n k]=get_dim(xTemp);
for i = 1:n
    %datTemp = data('tmp',[xTemp(:,1:(i-1)),xTemp(:,(i+1):size(xTemp,2))],[]);
   datTemp = get(xTemp,[],[1:(i-1) (i+1):n]);
   kerMai = get_kernel(supVec.child,datTemp,datTemp);      
    weights(i) = (alphaTemp'*kerMai*alphaTemp);
  i
end;   
weights=max(weights)-weights;
