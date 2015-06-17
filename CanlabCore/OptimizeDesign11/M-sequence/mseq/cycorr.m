function ccorr=cycorr(vect1,vect2)
% ccorr=cycorr(vect)
% calculates "cyclic" correlation

if nargin==1, vect2=vect1; end;

vect1=vect1(:);
vect2=vect2(:);

c=vect1';

n=length(c);
ccorr=zeros(1,n);
for i=1:n
   ccorr(i)=c*vect2;
   c=[c(2:end) c(1)];
end;
