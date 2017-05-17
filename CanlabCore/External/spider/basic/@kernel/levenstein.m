function K = levenstein(kern,dat1,dat2,ind1,ind2,kerParam)
Xf = get_x(dat1);
Xs = get_x(dat2);
X1 = {};
X2 = {};
for i = 1:size(Xf,1)
  tmp = Xf(i,:);
  X1{i} = tmp(tmp > 0);
end
for i = 1:size(Xs,1)
  tmp = Xs(i,:);
  X2{i} = tmp(tmp > 0);
end


%---calculating the distance matrix----
D= [];
for i = 1:length(ind1)
    disp(['D step ' num2str(i)])
    for j = 1:length(ind2)
        D(i,j) = d(X1{i},X2{j});
    end
end
%---get kernel matrix from distance matrix--
clear X1 X2 % free memory
K = [];
m = size(D);

const = prod(m.^(-1))*sum(sum(D.^2));
row_sum = m(1)^(-1)*sum(D.^2,1);
col_sum = m(2)^(-1)*sum(D.^2,2);
for i = ind1
    disp(['D->K step ' num2str(i)])
    for j = ind2
        K(i,j) = 0.5*(D(i,j)^2 - col_sum(i) - row_sum(j) + const);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = d(a,b)
    %----------initializing-------------
    C = zeros(length(a)+1,length(b)+1);
    C(:,1) = [0:length(a)]';
    C(1,:) = [0:length(b)];
    %------compute distance-------------
    for i = 2:size(C,1)
        for j = 2:size(C,2)
            delta = 1-abs(sign(a(i-1)-b(j-1)));
            C(i,j) = min([C(i-1,j)+1,C(i,j-1)+1,C(i-1,j-1) + 1 - delta]);  
        end
    end
    
    s = C(end,end);