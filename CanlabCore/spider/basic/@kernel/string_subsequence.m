function K = string_subsequence(kern,dat1,dat2,ind1,ind2,kerParam)
% --- conversion to sequences------------------
D= []; X1 = get_x(dat1); X2 = get_x(dat2);
% for i = 1:size(Xf,1), tmp = Xf(i,:);  X1{i} = tmp(tmp > 0); end
% for i = 1:size(Xs,1),  tmp = Xs(i,:);  X2{i} = tmp(tmp > 0);end

n = kerParam{1};
lambda = kerParam{2};
% disp(['Calculating subsequence kernel with lambda=' num2str(lambda) ' for all subsequences of length ' num2str(n)])

K = [];
for i = 1:length(X2)
%     disp(['Calculating line ' num2str(i)])
    for j = 1:length(X1)
       K(i,j) =  k(X1{j},X2{i},n,lambda);
    end
end
%------------------------------------------------------
function kret = k(s,t,n,lambda)
        if min(length(s),length(t)) < n
            kret = 0; 
        else
            kret = k(s(1:end-1),t,n,lambda);
            j = find(t == s(end));
            for i = 1:length(j)
                kret = kret + k_prime(s(1:end-1),t(1:j(i)-1),lambda,n-1)*lambda^2;   
            end
        end
        
%----------------------------------------------------------
function ret = k_prime(s,t,lambda,ind)
    if ind == 0
        ret = 1; return;
    elseif  min(length(s),length(t)) < ind
        ret = 0; return;
    else
        ret = lambda*k_prime(s(1:end-1),t,lambda,ind);
        j = find(t == s(end));
        for i = 1:length(j)
          ret = ret + k_prime(s(1:end-1),t(1:j(i)-1),lambda,ind-1)*lambda^(length(t)-j(i)+2);   
        end
    end

    