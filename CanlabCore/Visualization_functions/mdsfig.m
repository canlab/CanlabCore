function [f1,hh] = mdsfig(varargin)
% ::
%
%    [f1,hh] = mdsfig(pc,names,classes,linemat,cordir)
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
%    mdsfig(pc,clus,names,sigmat);
%
%    x = randn(5,2);
%    y = [1 1 1 2 2]';
%    c = eye(5); c(1,2) = 1; c(1,4) = -1; c(2,3) = 1;
%    mdsfig(x,y,[],c)
%
% :Example2: using data output from mvroi.m
% ::
%
%    DAT.coords = CLUSTER.Gs(:,1:2); DAT.names = CLUSTER.names; DAT.classes =
%    CLUSTER.classes; DAT.lines = DATA.CORRELS.AVGSTATS.sigmat; 
%    figure;
%    mdsfig(DAT.coords,DAT.names,DAT.classes,DAT.lines);

hold off;
f1 = [];
names = [];

% some defaults
if length(varargin) < 1
    error('must input at least group space');
end
X=varargin{1};
clus=ones(1,size(X,1));
for i = 1:size(X,1);
        names{i} = ['ROI' num2str(i)];
end;
linemat=zeros(length(clus),length(clus));
cordir=ones(length(clus),length(clus));
if length(varargin) > 1, names = varargin{2}; end
if length(varargin) > 2, clus = varargin{3}; end
if length(varargin) > 3, linemat = varargin{4}; end
if length(varargin) > 4, cordir = varargin{5}; end

%ignore zeros
X=X(sum(X,2)~=0,:);

if ~isempty(names)
    if length(names) == max(clus)
    else
        names={names{find(clus~=0)}};
    end
end

clus=clus(find(clus~=0));

%%% remove extra dimensions for plotting;
if size(X,2)>3;
    disp('removing last dimensions for plotting');
    X=X(:,1:3);
end

%if isempty(names);
%    for i = 1:size(X,1);
%        names{i} = ['R' num2str(i)];
%    end;
%end

if size(X,2)>2;
    disp('switching to 3D plotting');
    mdsfig_3d(X,names,clus,linemat,cordir);
    return
end

hold on;
%set(gca,'FontSize',10);
colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
while max(clus) > length(colors), colors = [colors colors];,end

[dummy,whfirst] = unique(clus); 
hh = [];

if length(names) == length(whfirst), 
    for j = 1:length(whfirst)
        hh(j) = plot(X(whfirst(j),1),X(whfirst(j),2),colors{clus(whfirst(j))},'MarkerSize',10,'MarkerFaceColor',colors{clus(whfirst(j))}(1));
    
    end
end
    
% Add lines

drawlines(X,linemat,0,{'Positive' 'Negative'});


% Plot points
for j = 1:size(X,1)
     h = plot(X(j,1),X(j,2),colors{clus(j)},'MarkerSize',10,'MarkerFaceColor',colors{clus(j)}(1));
    
     if any(j == whfirst), hh(end+1) = h;, end
        

end

if length(names) == length(whfirst), 
    legend(hh,names);,

else

    % Add text

    if ~isempty(names);
        for i = 1:size(X,1)
            text(X(i,1)+.025*[max(X(:,1))-min(X(:,1))],X(i,2)+.0,names{i},'Color','k','FontSize',12,'FontWeight','bold');

        end
    end
end


xlabel('Component 1'),ylabel('Component 2')
set(gca,'XTickLabel',{[]});  
set(gca,'YTickLabel',{[]});  
axis equal;

return





function [hhp,hhn] = drawlines(pc,sigmat,sigcol,legmat,varargin)
    
% sigcol is a scalar between zero and one, lower is 'more salient'
hold on

% linew
lw = 2;

color(1,:) = [0 0 0];   % first line color, 'positive'
color(2,:) = [0 .5 1];   % second line color, 'negative'
style{1} = '-';         % first line style
style{2} = '-';         % second line style

if length(varargin) > 0, color = varargin{1}; end
if length(varargin) > 1, style = varargin{2}; sigcol = 0; end

hhp=[];hhn=[];
for i=1:size(pc,1);
     for j=1:size(pc,1);
         if sigmat(i,j) > 0
             hhp = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
             set(hhp,'Color',color(1,:) + [1 1 1] * abs(sigcol),'LineStyle',style{1},'LineWidth',lw - (lw*sigcol)-1)
         elseif sigmat(i,j) < 0
             hhn = line([pc(i,1) pc(j,1)],[pc(i,2) pc(j,2)]);
             set(hhn,'Color',color(2,:) + abs([sigcol sigcol 0]),'LineStyle',style{2},'LineWidth',lw - (lw*sigcol)-1)
         end
     end
end
legend ([hhp hhn],legmat);
return
