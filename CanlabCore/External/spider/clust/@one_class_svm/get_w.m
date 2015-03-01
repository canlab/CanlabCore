
function w = get_w(sv) 
   
% w = get_w() 
% returns value of w from svm 
% AE: i have changed this procedure to include non linear feature
% selection with RFE. When non-linear, it gives  a score vector for each
% feature (used in RFE) which is not the margin but is sorted as the margin
if strcmp(sv.ker{1},'linear'),
  
  w=sv.alpha'*get_x(sv.Xsv);
  
else
    Xtmp = get_x(sv.Xsv);
    alphatmp = sv.alpha;
    
    if strcmp(sv.ker{1},'rbf'),
        
      % compute the kernel matrix for all components
      K = Xtmp*Xtmp';
      Kdn = sum(Xtmp.^2,2);
      Kn = sum(Xtmp.^2,2);
      K = ones(size(Xtmp,1),1)*Kdn' + Kn*ones(1,size(Xtmp,1)) - 2*K;
      K = exp(-K/(2*sv.ker{2}));
      % compute the margin when one component is removed
       for i = 1:size(Xtmp,2),
           Ki = Xtmp(:,i)*ones(1,size(Xtmp,1)) - ones(size(Xtmp,1),1)*Xtmp(:,i)';
           Ki = Ki.^2;
           Ki = Ki/(2*sv.ker{2}); 
           Ki = exp(Ki);
           w(i) = alphatmp'*(K.*Ki)*alphatmp;
       end;
             
    elseif strcmp(sv.ker{1},'poly'),
    
        % compute the margin when one component is removed
        Ktmp = Xtmp*Xtmp';
        for i = 1:size(Xtmp,2),
           Ki = Xtmp(:,i)*Xtmp(:,i)';
           K_i = (Ktmp - Ki+1).^(sv.ker{2});           
           w(i) = alphatmp'*K_i*alphatmp;
	end;
    end;% if strcmp(...,'rbf')
    w=max(w)-w;  %% make largest the smallest --- wrong way round!
    
end;% if strcmp(...,'linear')
