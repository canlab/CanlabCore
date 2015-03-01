function d=make_data(a)
X=[];
Y=[];

N=a.l;
name=a.dist;

% find the installation path
try 
	p=which(get_name(a));
	ii=findstr(p,['@',get_name(a)])+length(['@',get_name(a)]);
	img=double(imread( [p(1:ii-1),filesep,name,'.jpg']))/255;
catch
	img=double(imread( [name,'.jpg']))/255;
end

[dx,dy,c]=size(img);

%imshow(img)
n=1;
while( n<=N)
    u=ceil(rand(1,2).*[dx dy]);
    
    p= img(u(1),u(2),1)+img(u(1),u(2),3);
    if(p~=0)    
        if ( p>rand )
            X=[X; u];
            if(img(u(1),u(2),1)>img(u(1),u(2),3))
                Y=[Y;1];
            else
                Y=[Y;-1];
            end
            n=n+1;
        end
    end
    
    
end

%data is in unit square
d=data(X/128 -1,Y);