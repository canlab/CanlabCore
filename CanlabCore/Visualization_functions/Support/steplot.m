function steplot(x,plottype,inwid);
%function steplot(x,plottype,inwid);
% x must be a subs*data matrix. if 3d, will plot 2nd dim as separate lines.
% plottype can be 'b' (bar) or 'l' (line);
% inwid (optional) is the width of the se bars (default 1/10 of bar width).

if nargin<3
    inwid=10;
end

if nargin<2
    plottype='l';
end

if nargin>1;
    if plottype~='l' & plottype~='b'  & plottype~='p' & plottype~='x';
        disp('plottype must be line (l) or bar (b) or point or xy');
    end
end


if ndims(x)>3
    error('cant plot a 4d matrix');
end

if ndims(x)==2;

np=size(x,2);    
% get standard errors
steX=(std(x,[],1))/sqrt(size(x,1));     %standard error is std/root n
 
%mean
mx=squeeze(mean(x));

%reshape
%figure;
if plottype=='l';
plot(mx,'linewidth',4);
elseif plottype=='p';
plot(mx,'o','color','k','markerfacecolor','k','markersize',10);
elseif plottype=='b';    
bar(mx);
end
%   set(gca,'Ylim',[ 0 0.5]);

hold on;

for n=1:np;  
    line([n,n],[mx(n)-(steX(n)/2),mx(n)+(steX(n)/2)],'color','k','linewidth',2);
    line([n-1/inwid,n+1/inwid],[mx(n)-(steX(n)/2),mx(n)-(steX(n)/2)],'color','k','linewidth',2);
    line([n-1/inwid,n+1/inwid],[mx(n)+(steX(n)/2),mx(n)+(steX(n)/2)],'color','k','linewidth',2);
end

elseif ndims(x)==3;
    
if plottype=='l';    
bw=1;        %linewidth
np=size(x,3);
X=x(:,:);        
steX=std(X,[],1)/sqrt(size(x,1));
steX=reshape(steX,size(x,2),np);

mx=squeeze(mean(x));
%figure;
plot(mx','linewidth',3);
hold on;

for l=1:size(x,2);
for n=1:np;  
    line([n,n],[mx(l,n)-(steX(l,n)/2),mx(l,n)+(steX(l,n)/2)],'color','k','linewidth',bw);
    line([n-1/inwid,n+1/inwid],[mx(l,n)-(steX(l,n)/2),mx(l,n)-(steX(l,n)/2)],'color','k','linewidth',bw);
    line([n-1/inwid,n+1/inwid],[mx(l,n)+(steX(l,n)/2),mx(l,n)+(steX(l,n)/2)],'color','k','linewidth',bw);
end
end

elseif plottype=='b';
bw=2;
mx=squeeze(mean(x));
h=bar(mx');
h1=get(gca,'Children');a=get(h1);;
np=size(x,3);
X=x(:,:);        
steX=std(X,[],1)/sqrt(size(x,1));
steX=reshape(steX,size(x,2),np);

for l=1:size(x,2);
        a1=get(a(l).Children);a1.XData;
        wherex={a1.XData};wherey={a1.YData};
         whereX=wherex{1};whereY=wherey{1};
        wid=(whereX(3)-whereX(1))./inwid;
        wX=(whereX(3,:)+whereX(1,:))./2;
        wY=(whereY(2,:));
        for n=1:np;          
        line([wX(n) wX(n)],[wY(n)-steX(l,n) wY(n)+steX(l,n)],'color','k','linewidth',bw);
        line([wX(n)-wid wX(n)+wid],[wY(n)-steX(l,n) wY(n)-steX(l,n)],'color','k','linewidth',bw);
        line([wX(n)-wid wX(n)+wid],[wY(n)+steX(l,n) wY(n)+steX(l,n)],'color','k','linewidth',bw);        
     end
end


end
end
    