function [alpha,b,Xsv]= lsvm(optcell,X,Y),
% Interface to the libsvm 
%
% run the different algos of libsvm
%
%
disp('building the training set...');
fid = fopen('templibsvmxx.txt','r');
if fid==-1,
[m,d] = size(X);
fid = fopen('templibsvmxx.txt','w');
for i=1:m,
    fprintf(fid,'%f ',Y(i));
    for j=1:d,
        fprintf(fid,'%d:%f ',j,X(i,j)); 
    end;    
    fprintf(fid,'\n');
end;
else
    disp('WARNING: using the previous training set, if not erase the file templibsvmxx');
end; %%fid ==-1
fclose(fid);
m = size(X,1);
switch optcell.name
    
case 'C-SVC',
    disp('C-SVC');
    tmp = ' -s 0 ';
    
    kernel = optcell.kernel;
    
    C = optcell.C;
    
    C1 = optcell.C1;
    Cm1 = optcell.Cm1;
   
    
    if isempty(C1),
        C1 = 1;
        Cm1 = 1;
    end;
    
    
    C =min([10000,C]);
        
    tmp = [tmp , ' -c ' num2str(C)];
    tmp = [tmp , ' -w1 ' num2str(C1) ' -w-1 ' num2str(Cm1)];
    
    
    kernel = optcell.kernel;
    degree = optcell.deg;
    coef0 = optcell.coef0;
    
    if kernel<2,
        tmp = [tmp , ' -t ' num2str(kernel), ' -d ' num2str(degree) ' -r ' num2str(coef0)];
    else
        sigma = optcell.sigma;
        gam = 1/(2*sigma^2);
        tmp = [tmp, ' -t ' num2str(kernel) , ' -g ' num2str(gam)];
    end;    
        
case 'nu-SVC',
    disp('nu-SVC');
    tmp = ' -s 1 ';
    kernel = optcell.kernel;
    nu = optcell.nu;
    C1 = optcell.C1;
    Cm1 = optcell.Cm1;
    
    if isempty(C1),
        C1 = 1;
        Cm1 = 1;
    end;
    
           
    tmp = [tmp , ' -n ' num2str(nu)];
    tmp = [tmp , ' -w1 ' num2str(C1) ' -w-1 ' num2str(Cm1)];
     
    kernel = optcell.kernel;
    degree = optcell.deg;
    coef0 = optcell.coef0;
    
    if kernel<2,
        tmp = [tmp , ' -t ' num2str(kernel), ' -d ' num2str(degree) ' -r ' num2str(coef0)];
    else
        sigma = optcell.sigma;
        gam = 1/(2*sigma^2);
        tmp = [tmp, ' -t ' num2str(kernel) , ' -g ' num2str(gam)];
    end;
    
    
    
case 'eps-SVR',
    disp('eps-SVR');
     tmp = ' -s 3 ';
    kernel = optcell.kernel;
    C = optcell.C;
    C1 = optcell.C1;
    Cm1 = optcell.Cm1;
    epsilon = optcell.epsilon;
    
    if isempty(C1),
        C1 = 1;
        Cm1 = 1;
    end;
    
    C =min([10000,C]);
        
           
    tmp = [tmp , ' -c ' num2str(C)];
    %tmp = [tmp , ' -w1 ' num2str(C1) ' -w-1 ' num2str(Cm1)];
    tmp = [tmp , ' -p ' num2str(epsilon)];
    
    kernel = optcell.kernel;
    degree = optcell.deg;
    coef0 = optcell.coef0;
    
    if kernel<2,
        tmp = [tmp , ' -t ' num2str(kernel), ' -d ' num2str(degree) ' -r ' num2str(coef0)];
    else
        sigma = optcell.sigma;
        gam = 1/(2*sigma^2);
        tmp = [tmp, ' -t ' num2str(kernel) , ' -g ' num2str(gam)];
    end;
    
    
    
case 'one-class-SVM',
    
    disp('one class SVM');
     tmp = ' -s 2 ';
    kernel = optcell.kernel;
    nu = optcell.nu;
    
    tmp = [tmp , ' -n ' num2str(nu)];
    
    
    kernel = optcell.kernel;
    degree = optcell.deg;
    coef0 = optcell.coef0;
    
    if kernel<2,
        tmp = [tmp , ' -t ' num2str(kernel), ' -d ' num2str(degree) ' -r ' num2str(coef0)];
    else
        sigma = optcell.sigma;
        gam = 1/(2*sigma^2);
        tmp = [tmp, ' -t ' num2str(kernel) , ' -g ' num2str(gam)];
    end;
    
    
case 'nu-SVR',
    
    disp('nu-SVR');
    tmp = ' -s 4 ';
    C1 = optcell.C1;
    Cm1 = optcell.Cm1;
    nu = optcell.nu;
    C = optcell.C;    
    
    tmp = [tmp , ' -c ' num2str(C)];
    tmp = [tmp , ' -n ' num2str(nu)];
    tmp = [tmp , ' -w1 ' num2str(C1) ' -w-1 ' num2str(Cm1)];
    
    
    
    kernel = optcell.kernel;
    degree = optcell.deg;
    coef0 = optcell.coef0;
    
    if kernel<2,
        tmp = [tmp , ' -t ' num2str(kernel), ' -d ' num2str(degree) ' -r ' num2str(coef0)];
    else
        sigma = optcell.sigma;
        gam = 1/(2*sigma^2);
        tmp = [tmp, ' -t ' num2str(kernel) , ' -g ' num2str(gam)];
    end;
        
end;
if isunix,
      c_string = 'Optimization/libsvm/svm-train';
else
      c_string = 'Optimization/libsvm/svmtrain.exe';
end;
todo = [spider_path c_string tmp ' templibsvmxx.txt modelfinlibsvmxx'];
if isunix,
unix(todo);
else
dos(todo);
end;
disp('building the model...');
[m,d] = size(X);
fid = fopen('modelfinlibsvmxx','r');
finish=0;
while ~finish,
TMP = fgetl(fid);
t = findstr(TMP,'rho');
if ~isempty(t),
    b = str2num(TMP(t(1)+4:length(TMP)));
end;
t = findstr(TMP,'SV');
if ~isempty(t),
    finish=1;
end;
end;
%alpha = zeros(m,1);
count=0; 
Xsv=[];
alpha=[];
while feof(fid)==0,
    alphatmp = fscanf(fid,'%f ',1);
    tmp =fscanf(fid,'%*d:%f ',d);
    Xsv = [Xsv;tmp'];
    alpha = [alpha;alphatmp];
    %for i=1:m,
    %    if norm(tmp'-X(i,:))<10^(-3),
    %        alpha(i) = alphatmp;
    %        break;
    %    end;
    % end; 
end;
% because libsvm compute -b and not b
b=-b;
fclose(fid);
if isunix,
      unix('rm -f modelfinlibsvmxx');
      unix('rm -f templibsvmxx.txt');
else
      dos('del modelfinlibsvmxx');
      dos('del templibsvmxx.txt');
end;
