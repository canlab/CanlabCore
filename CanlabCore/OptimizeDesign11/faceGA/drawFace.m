function drawFace(p)
%
% head
% eyes
% nose
% mouth
% eyebrows
% hair
%
% 1 - head width    0 1
% 2 - head height   0 1
% 3-5 head color    0 1
%
% 6     y eye level
% 7     eye spacing
% 8     eye line length
% 9     left eye line width
% 10    right eye line width
%
% 11    nose to eye dist
% 12    nose length
% 13    nose line width
%
% 14    mouth width
% 15    smile scale factor
% 16    smile displacement
% 17    smile line width

cla

% ------------------------------------------
% * draw head
% ------------------------------------------
[X,Y,Z]=ellipsoid(0,0,0,p(1),p(2),1);
surf(X,Y,Z,'EdgeColor','None','FaceColor',[p(3) p(4) p(5)])
view(0,90); axis off; set(gcf,'Color','w');axis([-1 1 -1 1])
hold on

% ------------------------------------------
% * draw eyes
% ------------------------------------------
y = p(6) * p(2);
s = (p(7) * p(1)) ./ 2;
plot3([-p(8) * .2 - s -s],[y y],[2 2],'k','LineWidth',p(9)*8)
plot3([s (p(8)*.2+s)],[y y],[2 2],'k','LineWidth',p(10)*8)


% ------------------------------------------
% * draw nose
% ------------------------------------------
nosestart = p(2) - y - (p(11) * p(2));
noselen = p(12) * (p(2) - nosestart);
plot3([0 0],[nosestart nosestart-noselen],[2 2],'k','LineWidth',p(13)*8)

% ------------------------------------------
% * draw mouth
% ------------------------------------------
wid = (p(14) * p(1)) ./ 2;
x = -wid:.05:wid;
ys = p(15) * 2 * x .^ 2;
ys = ys - (p(16)*p(2));
plot3(x,ys,2*ones(size(x)),'k','LineWidth',p(17)*5)

drawnow

return
