function ret = kernel_loss(algo,dat)

%
%               function res = kernel_loss(algo,data)
%
% Loss derived from kernel matrix (supplied via param). The Kernel should
% be of type 'custom' or 'custom_fast' (i.e. all matrix entries should be
% stored).
%
% Note:
% The Kernel-Matrix represents the inner products in the 'loss' space
% between the examples.
% In Terms of this kernel we must compare output and input of a data
% object, so d.X (input) and d.Y (output) are indexes into its matrix.



if strcmp(algo.param.ker,'custom')
   kerPAram=algo.param.kerparam;
 else
   global K; 
   kerParam=K;
end

indToTake=dat.index; 
lss=0; 
j=1;
for i=indToTake
 Xrow=dat.X(j,:); 
 which=i;

 if length(Xrow)>1
   dist=0;
   for q=1:length(Xrow)
     dd=kerParam(Xrow(q),Xrow(q))+kerParam(which,which)-2*kerParam(Xrow(q),which);
     dist=dist+sqrt(dd);
   end
   lss=lss+dist/length(Xrow);
 else
   lss=lss+sqrt(kerParam(Xrow,Xrow)+kerParam(which,which)-2*kerParam(Xrow,which));
 end  
 j=j+1;
end
lss=lss/length(indToTake);

ret=data([get_name(dat) ' -> kernel_loss=' num2str(lss,4) ],[],lss);






