function bp=backprojection(c);
% function backp=backprojection(c);
% calculates backprojections from eigenvalues (e) and Group Space
%(c.GroupSpace) which are the outputs of cmdscale.

eigenvals=diag(c.eigenvalues(1:c.ndims),0);
A=c.GroupSpace(:,1:c.ndims)*sqrt(eigenvals);
for n=1:size(c.data,3);
    mdata=c.data(:,:,n);
    backproj(:,:,n)=A'*(mdata*inv(mdata'*mdata))';       
end
bp=mean(backproj,3);
