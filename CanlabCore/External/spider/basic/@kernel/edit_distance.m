function K = edit_distance(kern,dat1,dat2,ind1,ind2,kerParam)
%---calculating the distance matrix----
K= [];
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


for i = 1:length(ind1)
    disp(['D step ' num2str(i)])
    for j = 1:length(ind2)
        K(i,j) = max(length(X1{ind1(i)}),length(X2{ind2(j)})) - d(X1{ind1(i)},X2{ind2(j)});
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