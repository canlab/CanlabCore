function libsvm_save(f,d)

% libsvm_save(filename,data_object)

[l n k]=get_dim(d);
x=get_x(d);
y=get_y(d);
ff=fopen(f,'w');

for i=1:l
    
%    fprintf(ff,'%d ',sign((y(i)==1)-0.5));
    fprintf(ff,'%d ', y(i) );
    r=[1:0.5:n];
    r(2:2:n*2)=x(i,:);
    fprintf(ff,'%d:%f ',r);
    fprintf(ff,'\n');
    
end
fclose(ff);
