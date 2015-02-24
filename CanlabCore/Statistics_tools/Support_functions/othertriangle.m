function x=othertriangle(y);

if y(1,1)==1;
    warning('doesnt work for DIS matrices');
end

x=y*0;
if ndims(y)==2;
    x=y+rot90(flipud(triu(y)),3);
end

if ndims(y)==3;
    for n=1:size(y,3);
        y1=y(:,:,n);
        x(:,:,n)=y1+rot90(flipud(triu(y1)),3);
    end
end

if ndims(y)==4;
    for n=1:size(y,3);
        for n1=1:size(y,4);        
            y1=y(:,:,n,n1);
            x(:,:,n,n1)=y1+rot90(flipud(triu(y1)),3);
        end
    end
end