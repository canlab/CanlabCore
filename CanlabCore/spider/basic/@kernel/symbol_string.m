function K = symbol_string(kern,dat1,dat2,ind1,ind2,kerParam)
% --- conversion to sequences------------------
D= []; Xf = get_x(dat1); Xs = get_x(dat2); X1 = {}; X2 = {};
for i = 1:size(Xf,1), tmp = Xf(i,:);  X1{i} = tmp(tmp > 0); end
for i = 1:size(Xs,1),  tmp = Xs(i,:);  X2{i} = tmp(tmp > 0);end

n = kerParam{1};
lambda = kerParam{2};

for i = 1:size(X1,1)
    for j = i:size(X2,1)
       K(i,j) =  k(X1{i},X2{j},n,lambda)          
    end
end

%------------------------------------------------------
function kret = k(s,t,n,lambda)
        if min(length(s),length(t)) < n
            kret = 0; return
        else
            kret = k(s(1:end-1),t,n,lambda);
            ind = find(t == s(end));
            for i = 1:length(ind)
                kret = kret + k_prime();   
            end
        end
        
%----------------------------------------------------------
function ret = k_prime(s,t,lambda,i)
    if n == 0
        ret = 1; return;
    elseif min(length(s),length(t)) < i 
        ret = 0; return;
    end

    