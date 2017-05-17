function circplot(a,d)
%
% This function is ONLY for 2D data! It plots the code book and the
% given examples their radii (for dist and distd/ c.f. help kvq).
% It also colors the datapoint if labels are specified. Labels are 
% mandatory if the kvq object is in mode 'discriminative' or 'shared'.
% If you don't want to use label information for plotting set the mode
% to 'standard'.
%

	

	
x=d.X;
ax1=min(x); ax2=max(x);

if strcmp(a.mode,'discriminative') || strcmp(a.mode,'shared')
    y = get_y(d);
    h=plot(x(y==1,1),x(y==1,2),'g.','Markersize',7); hold on;
    h=plot(x(y==-1,1),x(y==-1,2),'m.','Markersize',7); hold on;
    ykeep = get_y(a.keep);
    h=plot(a.keep.X(ykeep==1,1),a.keep.X(ykeep==1,2),'go'); hold on;
    h=plot(a.keep.X(ykeep==-1,1),a.keep.X(ykeep==-1,2),'mo'); hold on;
else
	h=plot(x(:,1),x(:,2),'k.','Markersize',7); hold on;
    h=plot(a.keep.X(:,1),a.keep.X(:,2),'ko'); hold on;
end

% h=plot(x(:,1),x(:,2),'k.','Markersize',7); hold on;
% % set(h,'LineWidth',2,'MarkerSize',7);
% h=plot(a.keep.X(:,1),a.keep.X(:,2),'ko'); hold on;

granul=100;

R = a.dist;
t = [0:2*pi/granul:2*pi]; y = sin(t); x=cos(t);

sv = a.keep.X;
for i=1:size(sv,1)
    plot(R*x+sv(i,1),R*y+sv(i,2),'k');
end

if strcmp(a.mode,'discriminative') || strcmp(a.mode,'shared')
    Rd = a.distd;
    for i=1:size(sv,1)
        plot(Rd*x+sv(i,1),Rd*y+sv(i,2),'r');
    end
end


% set(h,'LineWidth',2,'MarkerSize',7);


hold off
