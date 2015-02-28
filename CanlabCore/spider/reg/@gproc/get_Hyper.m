function weig = get_w(supVec,method) 
%
%                w = get_w()
%
% Returns value of w from svm 
%
% AE: i have changed this procedure to include non linear feature
% selection with RFE. When non-linear, it gives  a score vector for each
% feature (used in RFE) which is not the margin but is sorted as the margin

kerTemp = sv.child.kerparam;
if strcmp(sv.child.ker,'linear')
  weig=sv.alpha'*get_x(sv.Xsv);
  return;
end

if strcmp(sv.child.ker,'weighted_linear')
  numEx=get_dim(sv.Xsv);
  weig=sv.alpha'* (get_x(sv.Xsv) .*  repmat(kerTemp,numEx,1));
  return;
end
if strcmp(sv.child.ker,'rbf')  | strcmp(sv.child.ker,'poly'),
    xTemp = get_x(sv.Xsv);
    alphaTemp = sv.alpha;
    
    if strcmp(sv.child.ker,'rbf'),
        
      % compute the kernel matrix for all components
      K = xTemp*xTemp';
      Kdn = sum(xTemp.^2,2);
      Kn = sum(xTemp.^2,2);
      K = ones(size(xTemp,1),1)*Kdn' + Kn*ones(1,size(xTemp,1)) - 2*K;
      K = exp(-K/(2*kerTemp));
      % compute the margin when one component is removed
       for i = 1:size(xTemp,2),
           Ki = xTemp(:,i)*ones(1,size(xTemp,1)) - ones(size(xTemp,1),1)*xTemp(:,i)';
           Ki = Ki.^2;
           Ki = Ki/(2*kerTemp); 
           Ki = exp(Ki);
           weig(i) = (alphaTemp'*(K.*Ki)*alphaTemp);
       end;
    elseif strcmp(sv.child.ker,'poly'),
        % compute the margin when one component is removed
        Ktmp = xTemp*xTemp';
        for i = 1:size(xTemp,2),
           Ki = xTemp(:,i)*xTemp(:,i)';
           K_i = (Ktmp - Ki+1).^(kerTemp);           
           weig(i) = (alphaTemp'*K_i*alphaTemp);
    	end;    
    end;% if strcmp(...,'rbf')
    weig=max(weig)-weig;  %% make largest the smallest --- wrong way round!
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
    [numEx,vDim,oDim]=get_dim(supVec.Xsv);    
    xTemp = get_x(supVec.Xsv);
    alphaTemp = supVec.alpha;
    weig = ones(length(supVec.child.param,vDim));
    
    for j=1:length(supVec.child.kerparam),  
        ktmp = supVec.child;
        kerp = ktmp.kerparam;
        for i = 1:length(kerp{j}),
            kerptmp=kerp;
            kerptmp{j} = [kerp{j}(1:i-1),kerp{j}(i+1:length(kerp{j}))];    
            ktmp.kerparam=kerptmp;
            K_i = get_kernel(ktmp,supVec.Xsv,supVec.Xsv);            
            weig(j,kerp{j}(i)) = (alphaTemp'*K_i*alphaTemp);
        end; 
    end;
    weig = max(max(weig))-weig;
    return;
end;
for i = 1:size(xTemp,2),
    dtmp = data('tmp',[xTemp(:,1:(i-1)),xTemp(:,(i+1):size(xTemp,2))],[]);
    Ki = get_kernel(supVec.child,dtmp,dtmp);      
    weig(i) = (alphaTemp'*Ki*alphaTemp);
end;   
weig=max(weig)-weig;