function f1 = nmdsfig1D(pc,clus,names,varargin)
% ::
%
%    f1 = nmdsfig(pc,clus,names)
%
% Create the plot with stimulus coordinates
%
% :Inputs:
%
%   **pc:**
%        is objects x dimensions
%
%   **clus:**
%        is a vector of object classes
%
%   **names:**
%        is cell array of names for rows of pc (objects), or empty ([])
%
%   **sigm:**
%        is optional matrix of 1, -1, and 0 entries
%        signifies which pairs to connect with lines
%        positive elements are solid lines, negative elements are dashed
%
% :Examples:
% ::
%
%    nmdsfig(pc,clus,names,sigmat);
%
%    x = randn(5,2);
%    y = [1 1 1 2 2]';
%    c = eye(5); c(1,2) = 1; c(1,4) = -1; c(2,3) = 1;
%    nmdsfig(x,y,[],c)


if size(pc,2)~=1; %%% remove extra dimensions for plotting;
error('only works with 1D configs. use nmdsfig() instead');
end


f1 = figure; set(gca,'FontSize',18),hold on
colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while max(clus) > length(colors), colors = [colors colors];,end


% Draw Lines
 
%if length(varargin) > 0; 
%    names  = varargin{1}; 
%if length(varargin) > 1
%    pmap  = varargin{2};
%else pmap=(sigmat*0)+1
end

if isempty(names);
    for i = 1:size(pc,1); 
        names{i} = ['V' num2str(i)];
    end;
end
%for i=1:size(pc,1);
%     for j=1:size(pc,1);
%         if sigmat(i,j) > 0
%             hh = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
%             set(hh,'Color',[.7 .7 .7],'LineWidth',fix(pmap(i,j)/10)+1)
%         elseif sigmat(i,j) < 0
%             hh = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
%             set(hh,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',fix(pmap(i,j)/10)+1)
%         end
%     end
%end

% Plot points
figure;
for j = 1:max(clus)
    hold on;
    plotpc=pc(clus == j);
     plot(plotpc,ones(1,length(plotpc)),colors{j},'MarkerFaceColor',colors{j}(1),'LineWidth',2,'MarkerSize' ,6)
end

for c=1:max(clus);
leg{c}=['clust',num2str(c)];
end
legend(leg);

xlabel('Component 1');
set(gca,'XTickLabel',{[]});  
return

